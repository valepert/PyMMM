from __future__ import division
from upwind import trasporto_upwind
from numpy import absolute, max, sqrt, sum, zeros


def errore_trasporto_upwind(n_iter, M, N, xmin, xmax, tinz, tfin, vel, ic):
    """Analisi errore trasporto upwind"""
    dx = zeros(n_iter)
    e1 = zeros(n_iter)
    e2 = zeros(n_iter)
    esup = zeros(n_iter)
    uexac = zeros([M + 1])

    dx[0] = (xmax - xmin) / M
    (u, x, t) = trasporto_upwind(M, N, xmin, xmax, tinz, tfin, vel, ic)
    uexac[0:M + 1] = ic(vel(x, t[-1]))

    e1[0] = dx[0] * sum(absolute(u[:, -1] - uexac))
    e2[0] = sqrt(dx[0] * sum(absolute(u[:, -1] - uexac) ** 2))
    esup[0] = max(abs(u[:, -1] - uexac))

    for i in range(1, n_iter):
        print i
        M *= 2
        dx[i] = (xmax - xmin) / M
        (u, x, t) = trasporto_upwind(M, N, xmin, xmax, tinz, tfin, vel, ic)
        uexac = zeros([M + 1])
        uexac[0:M + 1] = ic(vel(x, t[-1]))
        e1[i] = dx[i] * sum(absolute(u[:, -1] - uexac))
        e2[i] = sqrt(dx[i] * sum(absolute(u[:, -1] - uexac) ** 2))
        esup[i] = max(abs(u[:, -1] - uexac))

    return e1, e2, esup