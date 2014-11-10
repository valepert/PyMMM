from __future__ import division
from numpy import array, zeros


def trasporto_upwind(M, N, xmin, xmax, tinz, tfin, vel, ic):
    dx = (xmax - xmin) / M
    dt = (tfin - tinz) / N

    x = xmin + array(range(0, M + 1)) * dx
    t = array(range(0, N + 1)) * dt

    u = zeros([M + 1, N + 1])
    a = zeros([M + 1, N + 1])

    u[:, 0] = ic(x)

    for n in range(0, N + 1):
        a[0:M + 1, n] = vel(x, t[n])

    # CFL
    cfl = abs(a).max() * dt / dx

    if cfl > 1:
        print "CFL " + str(cfl)
        raise ArithmeticError
    else:
        for n in range(0, N):
            for j in range(1, M):
                nu = a[j, n] * dt / dx
                u[j, n + 1] = u[j, n] - 0.5 * nu * (u[j + 1, n] - u[j - 1, n]) + \
                              0.5 * abs(nu) * (u[j + 1, n] - 2 * u[j, n] + u[j - 1, n])

            u[0, n + 1] = 0.0
            u[M, N] = 0.0

    return u, x, t