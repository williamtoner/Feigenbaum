"""Low-precision prototype: degree-z Feigenbaum fixed point and delta via Chebyshev collocation.

g(x) = f(|x|^z), t = |x|^z in [0,1], lam = f(1) < 0, s = lam^z (z even)
F(t) = lam*f(t) - f(f(s t)^z) = 0;  alpha = -1/lam
delta = top eigenvalue of (L H)(t) = (1/lam) [ z u^(z-1) f'(u^z) H(s t) + H(u^z) ],  u = f(s t)
f(t) = sum' c_j T_j(2t-1)  (c_0 halved)
"""
import sys
import mpmath as mp


def nodes(n):
    return [(1 + mp.cos((2 * i - 1) * mp.pi / (2 * n))) / 2 for i in range(1, n + 1)]


def cheb_vals(n, t):
    """phi_k(t) = T_k(2t-1) for k < n, with phi_0 halved."""
    y = 2 * t - 1
    T = [mp.mpf(1), y]
    for _ in range(2, n):
        T.append(2 * y * T[-1] - T[-2])
    T = T[:n]
    T[0] = T[0] / 2
    return T


def cheb_dvals(n, t):
    """d/dt phi_k(t) = 2 T_k'(2t-1), via U_{k-1}: T_k' = k U_{k-1}."""
    y = 2 * t - 1
    U = [mp.mpf(1), 2 * y]
    for _ in range(2, n):
        U.append(2 * y * U[-1] - U[-2])
    return [mp.mpf(0)] + [2 * k * U[k - 1] for k in range(1, n)]


def ev(c, t):
    return mp.fsum(a * p for a, p in zip(c, cheb_vals(len(c), t)))


def evd(c, t):
    return mp.fsum(a * p for a, p in zip(c, cheb_dvals(len(c), t)))


def residual_and_jac(c, z, ts):
    n = len(c)
    lam = ev(c, 1)
    s = abs(lam) ** z
    dsdl = z * abs(lam) ** (z - 1) * mp.sign(lam)
    F, J = [], mp.matrix(n, n)
    for i, t in enumerate(ts):
        u = ev(c, s * t)
        w = u ** z
        fw = ev(c, w)
        fpw = evd(c, w)
        fpst = evd(c, s * t)
        ft = ev(c, t)
        F.append(lam * ft - fw)
        pt, pw, pst = cheb_vals(n, t), cheb_vals(n, w), cheb_vals(n, s * t)
        gu = fpw * z * u ** (z - 1)  # d f(u^z)/du
        for k in range(n):
            dl = mp.mpf(1) / 2 if k == 0 else mp.mpf(1)  # phi_k(1)
            du = pst[k] + fpst * t * dsdl * dl
            J[i, k] = dl * ft + lam * pt[k] - pw[k] - gu * du
    return F, J, lam


def solve_fixed_point(z, n, c0=None, iters=60, verbose=False):
    ts = nodes(n)
    c = list(c0) if c0 else [mp.mpf("0.6"), mp.mpf("-0.7")]
    c = (c + [mp.mpf(0)] * n)[:n]
    for it in range(iters):
        F, J, lam = residual_and_jac(c, z, ts)
        r = max(abs(x) for x in F)
        if verbose:
            print(f"  it {it} |F|={mp.nstr(r, 3)} alpha={mp.nstr(-1 / lam, 20)}")
        if r < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
        dc = mp.lu_solve(J, mp.matrix(F))
        c = [a - d for a, d in zip(c, dc)]
    return c, -1 / ev(c, 1)


def delta(c, z):
    n = len(c)
    ts = nodes(n)
    lam = ev(c, 1)
    s = abs(lam) ** z
    M, P = mp.matrix(n, n), mp.matrix(n, n)
    for i, t in enumerate(ts):
        u = ev(c, s * t)
        w = u ** z
        a = z * u ** (z - 1) * evd(c, w)
        pst, pw, pt = cheb_vals(n, s * t), cheb_vals(n, w), cheb_vals(n, t)
        for k in range(n):
            M[i, k] = (a * pst[k] + pw[k]) / lam
            P[i, k] = pt[k]
    E = mp.eig(mp.inverse(P) * M, left=False, right=False)
    return max(E, key=lambda e: abs(e))


if __name__ == "__main__":
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    c = None
    for z in [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8]:
        for n in [10, 20, 30, 40]:
            c, a = solve_fixed_point(z, n, c)
        if z in (2, 4, 6, 8):
            d = delta(c, z)
            print(f"z={z} n={n} alpha={mp.nstr(a, 30)} delta={mp.nstr(d, 30)}")
