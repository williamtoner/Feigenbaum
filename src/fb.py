"""High-precision Feigenbaum constants (alpha, delta) for degree-z unimodal maps, z even.

Formulation (t = |x|^z in [0,1], g(x) = f(t), f(t) = sum_k a_k T_k(2t-1)):
    lam = f(1) < 0,  s = lam^z,  u_i = f(s t_i),  w_i = u_i^z
    fixed point:  F_i = lam f(t_i) - f(w_i) = 0              (alpha = -1/lam)
    eigenpair:    G_i = (1/lam)[A_i H(s t_i) + H(w_i)] - delta H(t_i) = 0,
                  A_i = z u_i^(z-1) f'(w_i),  H = sum_k h_k T_k(2t-1),  h_0 = 1
Collocation at the n Chebyshev nodes of [0,1].

Full-precision arithmetic: Python ints as fixed point with P fractional bits (gmpy2 mpz).
Correction: x <- x - J^{-1} F, where J^{-1} F is solved by iterative refinement with a float64 LU
and an exact fixed-point Q-bit Jacobian (FLINT fmpz_mat), so each outer step gains ~Q bits.
Precision/node count ramp in stages.
"""
import argparse
import json
import math
import multiprocessing as mp_
import os
import sys
import time

import flint
import gmpy2
import mpmath as mp
import numpy as np
import scipy.linalg as sla
from gmpy2 import mpz

LOG2_10 = math.log2(10)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------- fixed point helpers

def fx_from_mpf(x, P):
    return mpz(mp.nint(mp.ldexp(x, P)))


def fmul(a, b, P):
    return (a * b) >> P


def fpow(a, z, P):
    r = a
    for _ in range(z - 1):
        r = (r * a) >> P
    return r


def clenshaw(a, y2, P):
    """sum_k a_k T_k(y), y2 = 2y in fixed point; a list of mpz."""
    b1 = mpz(0)
    b2 = mpz(0)
    for k in range(len(a) - 1, 0, -1):
        b1, b2 = a[k] + ((y2 * b1) >> P) - b2, b1
    return a[0] + ((y2 * b1) >> (P + 1)) - b2


def cheb_deriv(a, P):
    """Coefficients of df/dt where f(t) = sum a_k T_k(2t-1) (all with the 'full' T_0 convention)."""
    n = len(a)
    d = [mpz(0)] * (n + 1)
    for k in range(n - 1, 0, -1):
        d[k - 1] = d[k + 1] + 2 * k * a[k]
    d[0] = d[0] >> 1  # sum' convention -> full T_0 coefficient
    return [2 * x for x in d[:n]]  # chain rule dt: y = 2t-1


def nodes(n, P):
    """t_i = (1 + cos((2i-1) pi / 2n)) / 2 as fixed point ints."""
    with mp.workdps(int(P / LOG2_10) + 30):
        return [fx_from_mpf((1 + mp.cos((2 * i - 1) * mp.pi / (2 * n))) / 2, P) for i in range(1, n + 1)]


# ---------------------------------------------------------------- worker tasks

def _fp_chunk(args):
    """Fixed-point residual pieces for a chunk of nodes."""
    a, lam, s, ts, z, P = args
    one2 = mpz(2) << P
    out = []
    for t in ts:
        st = (s * t) >> P
        u = clenshaw(a, 4 * st - one2, P)
        w = fpow(u, z, P)
        fw = clenshaw(a, 4 * w - one2, P)
        ft = clenshaw(a, 4 * t - one2, P)
        out.append(((lam * ft) >> P) - fw)
    return out


def _aux_chunk(args):
    """u_i, w_i, s t_i, A_i = z u^(z-1) f'(w) for the eigen problem (full precision)."""
    a, da, s, ts, z, P = args
    one2 = mpz(2) << P
    out = []
    for t in ts:
        st = (s * t) >> P
        u = clenshaw(a, 4 * st - one2, P)
        w = fpow(u, z, P)
        A = (z * fpow(u, z - 1, P) * clenshaw(da, 4 * w - one2, P)) >> P if z > 1 else None
        out.append((st, w, A))
    return out


def _eig_chunk(args):
    h, ilam, delta, rows, P = args  # rows: (t, st, w, A)
    one2 = mpz(2) << P
    out = []
    for t, st, w, A in rows:
        hs = clenshaw(h, 4 * st - one2, P)
        hw = clenshaw(h, 4 * w - one2, P)
        ht = clenshaw(h, 4 * t - one2, P)
        out.append(((ilam * (((A * hs) >> P) + hw)) >> P) - ((delta * ht) >> P))
    return out


def _cheb_seq(y2, n, Q):
    """T_0..T_{n-1}(y) in Q-bit fixed point, y2 = 2y."""
    T = [mpz(1) << Q, y2 >> 1]
    for _ in range(2, n):
        T.append(((y2 * T[-1]) >> Q) - T[-2])
    return T[:n]


def _jac_rows(args):
    """Q-bit fixed-point Jacobian rows.
    mode 'fp': d F_i / d a_k
    mode 'eig': d G_i / d h_k (k>=1), column 0 replaced by d G_i / d delta = -H(t_i)."""
    mode, data, rows, n, z, Q = args
    one2 = mpz(2) << Q
    out = []
    if mode == "fp":
        lam, dsdl, a, da, s = data
        for t in rows:
            st = (s * t) >> Q
            u = clenshaw(a, 4 * st - one2, Q)
            w = fpow(u, z, Q)
            ft = clenshaw(a, 4 * t - one2, Q)
            g = (z * fpow(u, z - 1, Q) * clenshaw(da, 4 * w - one2, Q)) >> Q
            fpst_t_dsdl = (((clenshaw(da, 4 * st - one2, Q) * t) >> Q) * dsdl) >> Q
            Tt, Ts, Tw = _cheb_seq(4 * t - one2, n, Q), _cheb_seq(4 * st - one2, n, Q), _cheb_seq(4 * w - one2, n, Q)
            c0 = ft - ((g * fpst_t_dsdl) >> Q)  # k-independent part (d lam / d a_k = 1)
            out.append([int(c0 + ((lam * Tt[k]) >> Q) - Tw[k] - ((g * Ts[k]) >> Q)) for k in range(n)])
    else:
        ilam, delta, h = data
        for t, st, w, A in rows:
            Tt, Ts, Tw = _cheb_seq(4 * t - one2, n, Q), _cheb_seq(4 * st - one2, n, Q), _cheb_seq(4 * w - one2, n, Q)
            row = [((ilam * (((A * Ts[k]) >> Q) + Tw[k])) >> Q) - ((delta * Tt[k]) >> Q) for k in range(n)]
            row[0] = -clenshaw(h, 4 * t - one2, Q)
            out.append([int(x) for x in row])
    return out


# ---------------------------------------------------------------- solver

class Solver:
    def __init__(self, z, workers):
        self.z = z
        self.pool = mp_.get_context("fork").Pool(workers)
        self.W = workers
        flint.ctx.threads = workers

    def chunks(self, xs):
        k = max(1, math.ceil(len(xs) / (4 * self.W)))
        return [xs[i:i + k] for i in range(0, len(xs), k)]

    def pmap(self, fn, make_args, xs):
        res = self.pool.map(fn, [make_args(c) for c in self.chunks(xs)])
        return [y for r in res for y in r]

    # --- linear correction machinery
    def build_precond(self, rows_fn, n, Q):
        Jrows = rows_fn()
        Jint = flint.fmpz_mat(Jrows)
        Jf = np.array([[float(x) for x in r] for r in Jrows]) / 2.0 ** Q
        lu = sla.lu_factor(Jf)
        cond = np.linalg.cond(Jf) if n <= 3000 else float("nan")
        return Jint, lu, cond

    def lin_solve(self, Jint, lu, F, P, Q, inner_tol_bits):
        """Return Delta (P-bit fixed point) ~ J^{-1} F to ~inner_tol_bits relative accuracy.
        Iterative refinement: float64 LU for corrections, exact Q-bit integer Jacobian for residuals."""
        m = max(abs(x) for x in F)
        if m == 0:
            return [mpz(0)] * len(F)
        R = 64
        sh = m.bit_length() - Q  # Fq ~ 2^Q
        Fq = [x >> sh if sh >= 0 else x << -sh for x in F]
        target = [x << (Q + R) for x in Fq]  # J_int d_int = target  <=>  d_int = J^{-1} Fq 2^R
        tbits = max(abs(x) for x in target).bit_length()
        d = [mpz(0)] * len(F)
        r = target
        for it in range(500):
            rm = max(abs(x) for x in r)
            if rm == 0 or rm.bit_length() < tbits - inner_tol_bits:
                break
            e = rm.bit_length() - 60
            rf = np.array([float(x >> e) if e >= 0 else float(x << -e) for x in r])
            dd = sla.lu_solve(lu, rf)  # J^{-1} r 2^-e  (J in real units);  delta d_int = J^{-1} r / 2^Q
            k = e - Q
            inc = [mpz(int(v)) for v in np.rint(dd * 2.0 ** 20)]
            k -= 20
            d = [di + (x << k if k >= 0 else x >> -k) for di, x in zip(d, inc)]
            Jd = Jint * flint.fmpz_mat([[int(x)] for x in d])
            r = [t - mpz(int(Jd[i, 0])) for i, t in enumerate(target)]
        else:
            raise RuntimeError("inner refinement did not converge")
        self.inner_its = it
        k = sh - R
        return [x << k if k >= 0 else x >> -k for x in d]

    # --- fixed point
    def fp_residual(self, a, ts, P):
        lam = sum(a)
        s = fpow(lam, self.z, P)
        return self.pmap(_fp_chunk, lambda c: (a, lam, s, c, self.z, P), ts), lam

    def fp_jac(self, a, ts, n, Q, P):
        aq = [x >> (P - Q) for x in a]
        lam = sum(aq)
        s = fpow(lam, self.z, Q)
        dsdl = self.z * fpow(lam, self.z - 1, Q)
        da = cheb_deriv(aq, Q)
        tq = [t >> (P - Q) for t in ts]
        return lambda: self.pmap(_jac_rows, lambda c: ("fp", (lam, dsdl, aq, da, s), c, n, self.z, Q), tq)

    def solve_fp(self, a, n, P, Q, tol_bits, max_outer=10000):
        ts = nodes(n, P)
        t1 = time.time()
        Jint, lu, cond = self.build_precond(self.fp_jac(a, ts, n, Q, P), n, Q)
        log(f"  precond built in {time.time() - t1:.1f}s")
        log(f"  fp precond n={n} Q={Q} cond={cond:.2e}")
        for it in range(max_outer):
            t1 = time.time()
            F, lam = self.fp_residual(a, ts, P)
            t_res = time.time() - t1
            m = max(abs(x) for x in F)
            res_bits = P - m.bit_length() if m else P
            log(f"  fp it={it} n={n} P={P} |F|~2^-{res_bits} inner={getattr(self, 'inner_its', 0)} t_res={t_res:.1f}s")
            if res_bits >= tol_bits:
                return a, lam
            D = self.lin_solve(Jint, lu, F, P, Q, Q - 30)
            a = [x - y for x, y in zip(a, D)]
        raise RuntimeError("fp did not converge")

    # --- eigenproblem
    def eig_aux(self, a, ts, P):
        lam = sum(a)
        s = fpow(lam, self.z, P)
        da = cheb_deriv(a, P)
        aux = self.pmap(_aux_chunk, lambda c: (a, da, s, c, self.z, P), ts)
        ilam = (mpz(1) << (2 * P)) // lam
        return [(t, st, w, A) for t, (st, w, A) in zip(ts, aux)], ilam

    def eig_residual(self, h, ilam, delta, rows, P):
        return self.pmap(_eig_chunk, lambda c: (h, ilam, delta, c, P), rows)

    def solve_eig(self, a, h, delta, n, P, Q, tol_bits, max_outer=10000):
        ts = nodes(n, P)
        rows, ilam = self.eig_aux(a, ts, P)
        sh = P - Q
        rq = [(t >> sh, st >> sh, w >> sh, A >> sh) for t, st, w, A in rows]
        hq = [x >> sh for x in h]
        Jint, lu, cond = self.build_precond(
            lambda: self.pmap(_jac_rows, lambda c: ("eig", (ilam >> sh, delta >> sh, hq), c, n, self.z, Q), rq), n, Q)
        log(f"  eig precond n={n} Q={Q} cond={cond:.2e}")
        for it in range(max_outer):
            t1 = time.time()
            G = self.eig_residual(h, ilam, delta, rows, P)
            t_res = time.time() - t1
            m = max(abs(x) for x in G)
            res_bits = P - m.bit_length() if m else P
            log(f"  eig it={it} n={n} |G|~2^-{res_bits} inner={getattr(self, 'inner_its', 0)} t_res={t_res:.1f}s")
            if res_bits >= tol_bits:
                return h, delta
            D = self.lin_solve(Jint, lu, G, P, Q, Q - 30)
            delta = delta - D[0]
            h = [h[0]] + [x - y for x, y in zip(h[1:], D[1:])]
        raise RuntimeError("eig did not converge")


def rescale(v, P_old, P_new, n_new):
    sh = P_new - P_old
    v = [x << sh if sh >= 0 else x >> -sh for x in v]
    return (v + [mpz(0)] * n_new)[:n_new]


def bootstrap(z):
    """Low-precision start from the mpmath prototype (continuation in z)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import proto
    mp.mp.dps = 40
    c = None
    zs = [2] + [x / 2 for x in range(5, 2 * z + 1)]
    for zz in zs:
        for n in [10, 20, 30]:
            c, _ = proto.solve_fixed_point(zz, n, c)
    c, _ = proto.solve_fixed_point(z, 40, c)
    return [mp.mpf(x) for x in c]


def digits_of(xfx, P, nd):
    with mp.workdps(nd + 20):
        return mp.nstr(mp.ldexp(mp.mpf(int(xfx)), -P), nd + 5, strip_zeros=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=int, default=2)
    ap.add_argument("--digits", type=int, default=1100)
    ap.add_argument("--rate", type=float, default=None, help="digits per node (default by z)")
    ap.add_argument("--nfactor", type=float, default=1.0, help="multiply node counts (check runs)")
    ap.add_argument("--Q", type=int, default=400)
    ap.add_argument("--guard", type=int, default=96)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--no-delta", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    z = args.z
    rate = args.rate or {2: 1.55, 4: 0.95, 6: 0.78, 8: 0.66}.get(z, 0.6)
    S = Solver(z, args.workers)
    log(f"z={z} digits={args.digits} rate={rate} nfactor={args.nfactor} workers={args.workers}")
    t0 = time.time()

    # ---- bootstrap
    c = bootstrap(z)
    P0 = 160
    a = [fx_from_mpf(x, P0) for x in c]
    a[0] = a[0] >> 1  # sum' -> full convention
    P_prev = P0

    # stages: digits doubling
    stages = []
    d = 60
    while d < args.digits:
        d = min(2 * d, args.digits)
        stages.append(d)
    if not stages:
        stages = [args.digits]
    for D in stages:
        n = int(math.ceil((D + 30) / rate * args.nfactor)) + 10
        P = int(math.ceil(D * LOG2_10)) + args.guard
        a = rescale(a, P_prev, P, n)
        P_prev = P
        Q = min(args.Q, P)
        log(f"stage D={D} n={n} P={P}")
        a, lam = S.solve_fp(a, n, P, Q, tol_bits=P - args.guard // 2)
    lam_fx = sum(a)
    with mp.workdps(int(P / LOG2_10) + 20):
        alpha = -1 / mp.ldexp(mp.mpf(int(lam_fx)), -P)
        alpha_s = mp.nstr(alpha, args.digits + 20, strip_zeros=False)
    log(f"alpha done ({time.time() - t0:.0f}s): {alpha_s[:60]}...")
    result = {"z": z, "digits_target": args.digits, "n": n, "P": P, "nfactor": args.nfactor,
              "alpha": alpha_s, "t_alpha": time.time() - t0}
    if args.out:
        json.dump(result, open(args.out, "w"))

    if not args.no_delta:
        # bootstrap eigenpair at low precision with the mpmath prototype matrix
        import proto
        mp.mp.dps = 40
        nb = 30
        cb = [mp.ldexp(mp.mpf(int(x)), -P) for x in a[:nb]]
        cb[0] *= 2
        dl = proto.delta(cb, z).real
        # eigenvector: solve (M - dl P) v = 0 with v0 = 1 at low precision, via prototype matrices
        hb, dlt = eig_bootstrap(cb, z, dl)
        Pp = 160
        h = [fx_from_mpf(x, Pp) for x in hb]
        h = [mpz(1) << Pp] + h[1:]
        delta = fx_from_mpf(dlt, Pp)
        P_prev = Pp
        a_full, P_full = a, P
        for D in stages:
            n = int(math.ceil((D + 30) / rate * args.nfactor)) + 10
            P = int(math.ceil(D * LOG2_10)) + args.guard
            h = rescale(h, P_prev, P, n)
            delta = delta << (P - P_prev)
            P_prev = P
            af = rescale(a_full, P_full, P, len(a_full))[:n]
            Q = min(args.Q, P)
            log(f"eig stage D={D} n={n} P={P}")
            h, delta = S.solve_eig(af, h, delta, n, P, Q, tol_bits=P - args.guard // 2)
        with mp.workdps(int(P / LOG2_10) + 20):
            delta_s = mp.nstr(mp.ldexp(mp.mpf(int(delta)), -P), args.digits + 20, strip_zeros=False)
        log(f"delta done ({time.time() - t0:.0f}s): {delta_s[:60]}...")
        result.update({"delta": delta_s, "t_total": time.time() - t0})
        if args.out:
            json.dump(result, open(args.out, "w"))


def eig_bootstrap(c, z, dl):
    """Eigenvector (full-T0 convention, h0 = 1) of the collocation operator at low precision."""
    import proto
    n = len(c)
    ts = proto.nodes(n)
    lam = proto.ev(c, 1)
    s = abs(lam) ** z
    M = mp.matrix(n, n)
    for i, t in enumerate(ts):
        u = proto.ev(c, s * t)
        w = u ** z
        A = z * u ** (z - 1) * proto.evd(c, w)
        pst, pw, pt = proto.cheb_vals(n, s * t), proto.cheb_vals(n, w), proto.cheb_vals(n, t)
        for k in range(n):
            M[i, k] = (A * pst[k] + pw[k]) / lam - dl * pt[k]
    # fix v0 = 1 (halved-T0 basis in proto -> full basis: coefficient of T_0 is v0/2)
    Ms = mp.matrix(n, n - 1)
    rhs = mp.matrix(n, 1)
    for i in range(n):
        rhs[i] = -M[i, 0] * 2
        for k in range(1, n):
            Ms[i, k - 1] = M[i, k]
    v = mp.qr_solve(Ms, rhs)[0]
    return [mp.mpf(1)] + [v[k] for k in range(n - 1)], dl


if __name__ == "__main__":
    main()
