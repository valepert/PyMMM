from __future__ import division
from numpy import array, zeros


def trasporto_lw_bcp(M, N, xmin, xmax, tinz, tfin, vel, ic):
    dx = (xmax - xmin) / M
    dt = (tfin - tinz) / N

    x = xmin + array(range(0, M + 1)) * dx
    t = array(range(0, N + 1)) * dt

    u = zeros([M + 1, N + 1])
    u[:, 0] = ic(x)
    a = vel

    # CFL
    cfl = a * dt / dx
    # numero di Courant
    nu = a * dt / dx

    if cfl > 1:
        print "CFL " + str(cfl)
        raise ArithmeticError
    else:
        for n in range(0, N):
            u[0, n + 1] = 0.5 * nu * (nu + 1) * u[M - 1, n] + (1 - nu * nu) * u[0, n] - 0.5 * nu * (1 - nu) * u[1, n]

            for j in range(1, M):
                u[j, n + 1] = 0.5 * nu * (nu + 1) * u[j - 1, n] + (1 - nu * nu) * u[j, n] - 0.5 * nu * (1 - nu) * u[
                    j + 1, n]

            u[M, n + 1] = u[0, n + 1]

    return u, x, t