"""Independent second implementation (shares no code with fb.py / proto.py).

Differences from fb.py, on purpose:
  - basis: shifted Legendre polynomials P_k(2t-1)            (fb: Chebyshev T_k)
  - nodes: Chebyshev extrema (Gauss-Lobatto) on [0,1]         (fb: Chebyshev roots)
  - normalisation: explicit equation f(0) = 1                 (fb: implicit)
  - equation form: f(t) - f(w)/lam = 0                        (fb: lam f(t) - f(w) = 0)
  - eigenvector normalisation: Legendre coefficient h_0 = 1   (fb: Chebyshev c_0 = 1)
  - arithmetic: FLINT/Arb ball arithmetic (python-flint arb)  (fb: gmpy2 fixed point)
  - solver: full Newton with exact Jacobian, dense solve at working precision, precision doubling
            (fb: fixed float64/Q-bit preconditioner, linear convergence)
  - bootstrap: own continuation in z from f = 1 - 1.5 t, own small-n eigen solve via acb_mat.eig

Setting: g(x) = f(t), t = |x|^z (z even), lam = f(1), s = lam^z, u = f(s t), w = u^z.
  alpha = -1/lam;   delta = dominant eigenvalue of H -> (1/lam)[z u^(z-1) f'(w) H(s t) + H(w)].
"""
import argparse
import json
import math
import sys
import time

import flint
from flint import arb, arb_mat, acb_mat


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def legendre_all(y, n, deriv=False):
    """P_0..P_{n-1}(y) (and derivatives wrt y if deriv)."""
    P = [arb(1), y]
    for k in range(1, n - 1):
        P.append(((2 * k + 1) * y * P[k] - k * P[k - 1]) / (k + 1))
    P = P[:n]
    if not deriv:
        return P
    D = [arb(0), arb(1)]
    for k in range(1, n - 1):
        D.append(D[k - 1] + (2 * k + 1) * P[k])
    return P, D[:n]


def ev(a, x):
    """f(x) = sum a_k P_k(2x-1) and f'(x)."""
    P, D = legendre_all(2 * x - 1, len(a), deriv=True)
    f = arb(0)
    fp = arb(0)
    for ak, p, d in zip(a, P, D):
        f += ak * p
        fp += ak * d
    return f, 2 * fp


def lobatto(m):
    """m Chebyshev-extrema nodes on [0,1] including both endpoints: (1 - cos(pi i/(m-1)))/2."""
    return [(1 - arb.cos_pi(arb(i) / (m - 1))) / 2 for i in range(m)]


# ------------------------------------------------------------------ fixed point

def fp_system(a, z, nodes):
    """Residual vector and Jacobian. Unknowns a_0..a_{n-1}.
    Equations: f(0) - 1 = 0 and f(t_i) - f(w_i)/lam = 0 at n-1 nodes in (0,1]."""
    n = len(a)
    lam, _ = ev(a, arb(1))
    s = abs(lam) ** z  # lam < 0; z may be non-integer during the bootstrap continuation
    dsdl = -z * abs(lam) ** (z - 1)
    R = []
    J = []
    # normalisation row
    f0, _ = ev(a, arb(0))
    R.append(f0 - 1)
    J.append([arb(1) if k % 2 == 0 else arb(-1) for k in range(n)])
    for t in nodes:
        st = s * t
        u, fpst = ev(a, st)
        w = u ** z
        fw, fpw = ev(a, w)
        ft, _ = ev(a, t)
        R.append(ft - fw / lam)
        Pt = legendre_all(2 * t - 1, n)
        Ps = legendre_all(2 * st - 1, n)
        Pw = legendre_all(2 * w - 1, n)
        gu = fpw * z * u ** (z - 1)
        common = fw / lam ** 2 - gu * fpst * t * dsdl / lam  # k-independent part (d lam/d a_k = 1)
        J.append([Pt[k] - (Pw[k] + gu * Ps[k]) / lam + common for k in range(n)])
    return R, J


def newton_fp(a, z, n, target_digits, max_steps=12):
    a = [x.mid() for x in (list(a) + [arb(0)] * n)[:n]]
    nodes = lobatto(n)[1:]  # exclude t=0 (handled by the normalisation row)
    tol = arb(10) ** (-target_digits - 5)
    for step in range(max_steps):
        t0 = time.time()
        R, J = fp_system(a, z, nodes)
        rmax = max(abs(r.mid()) for r in R)
        log(f"  fp n={n} step={step} |R|={rmax.str(3, radius=False)} build={time.time() - t0:.0f}s")
        if rmax < tol:
            return a
        t0 = time.time()
        d = arb_mat(J).solve(arb_mat([[r] for r in R]), algorithm="approx")
        a = [(ak - d[k, 0]).mid() for k, ak in enumerate(a)]
        log(f"    solve {time.time() - t0:.0f}s")
    raise RuntimeError("fp Newton did not converge")


# ------------------------------------------------------------------ eigenproblem

def eig_aux(a, z, nodes):
    lam, _ = ev(a, arb(1))
    s = abs(lam) ** z
    rows = []
    for t in nodes:
        st = s * t
        u, _ = ev(a, st)
        w = u ** z
        _, fpw = ev(a, w)
        rows.append((t, st, w, z * u ** (z - 1) * fpw))
    return lam, rows


def eig_system(h, delta, lam, rows):
    """Unknowns h_1..h_{n-1}, delta (h_0 = 1). Equations at all n Lobatto nodes."""
    n = len(h)
    R, J = [], []
    for t, st, w, A in rows:
        Pt = legendre_all(2 * t - 1, n)
        Ps = legendre_all(2 * st - 1, n)
        Pw = legendre_all(2 * w - 1, n)
        Ht = sum((hk * p for hk, p in zip(h, Pt)), arb(0))
        Hs = sum((hk * p for hk, p in zip(h, Ps)), arb(0))
        Hw = sum((hk * p for hk, p in zip(h, Pw)), arb(0))
        R.append((A * Hs + Hw) / lam - delta * Ht)
        row = [(A * Ps[k] + Pw[k]) / lam - delta * Pt[k] for k in range(n)]
        row[0] = -Ht  # column 0 <-> delta
        J.append(row)
    return R, J


def newton_eig(a, h, delta, z, n, target_digits, max_steps=12):
    h = [x.mid() for x in (list(h) + [arb(0)] * n)[:n]]
    h[0] = arb(1)
    delta = delta.mid()
    a = [x.mid() for x in a]
    nodes = lobatto(n)
    lam, rows = eig_aux(a, z, nodes)
    tol = arb(10) ** (-target_digits - 5)
    for step in range(max_steps):
        t0 = time.time()
        R, J = eig_system(h, delta, lam, rows)
        rmax = max(abs(r.mid()) for r in R)
        log(f"  eig n={n} step={step} |R|={rmax.str(3, radius=False)} build={time.time() - t0:.0f}s")
        if rmax < tol:
            return h, delta
        t0 = time.time()
        d = arb_mat(J).solve(arb_mat([[r] for r in R]), algorithm="approx")
        delta = (delta - d[0, 0]).mid()
        h = [h[0]] + [(hk - d[k, 0]).mid() for k, hk in enumerate(h) if k > 0]
        log(f"    solve {time.time() - t0:.0f}s")
    raise RuntimeError("eig Newton did not converge")


def eig_bootstrap(a, z, n):
    """Dominant eigenpair of the collocated operator at small n via acb_mat.eig."""
    nodes = lobatto(n)
    lam, rows = eig_aux(a[:n] if len(a) >= n else a + [arb(0)] * (n - len(a)), z, nodes)
    M, Phi = [], []
    for t, st, w, A in rows:
        Pt = legendre_all(2 * t - 1, n)
        Ps = legendre_all(2 * st - 1, n)
        Pw = legendre_all(2 * w - 1, n)
        M.append([(A * Ps[k] + Pw[k]) / lam for k in range(n)])
        Phi.append(Pt)
    K = arb_mat(Phi).solve(arb_mat(M), algorithm="approx")  # coefficient-space operator
    ev_ = acb_mat(K).eig(algorithm="approx")
    best = max(ev_, key=lambda e: abs(e.real) if abs(e.imag.mid()) < 1e-10 else -1)
    d0 = best.real
    # eigenvector with h_0 = 1: solve (K - d0 I) h = 0 for h_1..h_{n-1}
    A = [[K[i, k] - (d0 if i == k else 0) for k in range(1, n)] for i in range(n)]
    rhs = [[-(K[i, 0] - (d0 if i == 0 else 0))] for i in range(n)]
    At = arb_mat(A).transpose()
    sol = (At * arb_mat(A)).solve(At * arb_mat(rhs), algorithm="approx")
    return [arb(1)] + [sol[k, 0] for k in range(n - 1)], d0


# ------------------------------------------------------------------ driver

def bootstrap_fp(z):
    flint.ctx.prec = 200
    a = [arb(0.25), arb(-0.75)]  # f = 1/4 - 3/4 (2t-1) = 1 - 1.5 t
    zs = [2.0] + [x / 4 for x in range(9, int(4 * z) + 1)]
    for zz in zs:
        for n in (6, 12, 20):
            a = newton_fp(a, zz, n, 40)
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=int, required=True)
    ap.add_argument("--digits", type=int, required=True)
    ap.add_argument("--rate", type=float, default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    flint.ctx.threads = args.threads
    z = args.z
    rate = args.rate or {2: 1.5, 4: 0.9, 6: 0.74, 8: 0.62}[z]
    t_start = time.time()
    log(f"alt z={z} digits={args.digits}")

    a = bootstrap_fp(z)
    stages = []
    d = 40
    while d < args.digits:
        d = min(int(d * 1.8), args.digits)
        stages.append(d)
    for D in stages:
        n = int(math.ceil((D + 25) / rate)) + 10
        flint.ctx.prec = int(D * 3.33) + 128
        log(f"fp stage D={D} n={n} prec={flint.ctx.prec}")
        a = newton_fp(a, z, n, D)
    lam, _ = ev(a, arb(1))
    alpha = -1 / lam
    res = {"z": z, "digits": args.digits, "n": n, "alpha": alpha.str(args.digits + 15, radius=False)}
    log(f"alpha {res['alpha'][:50]}... rad={alpha.rad()} ({time.time() - t_start:.0f}s)")
    json.dump(res, open(args.out, "w"))

    flint.ctx.prec = 200
    h, delta = eig_bootstrap([x for x in a], z, 24)
    log(f"eig bootstrap delta ~ {delta.str(15, radius=False)}")
    d = 40
    for D in stages:
        n = int(math.ceil((D + 25) / rate)) + 10
        flint.ctx.prec = int(D * 3.33) + 128
        log(f"eig stage D={D} n={n}")
        h, delta = newton_eig(a, h, delta, z, n, D)
    res["delta"] = delta.str(args.digits + 15, radius=False)
    res["t_total"] = time.time() - t_start
    log(f"delta {res['delta'][:50]}... ({res['t_total']:.0f}s)")
    json.dump(res, open(args.out, "w"))


if __name__ == "__main__":
    main()
