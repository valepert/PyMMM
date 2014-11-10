from __future__ import division
from numpy import array, diagflat, eye, ones, zeros
from numpy.linalg import solve


def trasporto_ei_bcp(M, N, xmin, xmax, tinz, tfin, vel, ic):
    dx = (xmax - xmin) / M
    dt = (tfin - tinz) / N

    x = xmin + array(range(0, M + 1)) * dx
    t = array(range(0, N + 1)) * dt

    u = zeros([M + 1, N + 1])
    u[:, 0] = ic(x)

    uold = zeros([M + 1])
    uold[0:M + 1] = u[:, 0]

    a = vel
    # CFL
    cfl = a * dt / dx
    # numero di Courant
    nu = a * dt / dx

    A = eye(M + 1) + 0.5 * nu * diagflat(ones([1, M]), 1) - 0.5 * nu * diagflat(ones([1, M]), -1)
    A[0, M] = -0.5 * nu
    A[M, 1] = 0.5 * nu

    if cfl > 1:
        print "CFL " + str(cfl)
        raise ArithmeticError

    else:
        for n in range(0, N + 1):
            unew = solve(A, uold.transpose()).transpose()
            uold = unew
            u[0:M + 1, n] = unew

    return u, x, t


def burgers_ei_bcp(M, N, xmin, xmax, tinz, tfin, ic):
    dx = (xmax - xmin) / M
    dt = (tfin - tinz) / N

    x = xmin + array(range(0, M + 1)) * dx
    t = array(range(0, N + 1)) * dt

    u = zeros([M + 1, N + 1])
    u[:, 0] = ic(x)

    uold = zeros([M + 1])
    uold[0:M + 1] = u[:, 0]

    for n in range(0, N + 1):
        # numero di Courant
        nu = max(abs(uold)) * dt / dx
        # CFL
        cfl = nu * dt / dx

        if cfl > 1:
            print "CFL " + str(cfl)
            raise ArithmeticError

        else:
            A = eye(M + 1) + 0.5 * dt * diagflat(uold[0:M], 1) / dx - 0.5 * dt * diagflat(uold[1:M + 1], -1) / dx
            A[0, M] = -0.5 * dt * uold[0] / dx
            A[M, 1] = 0.5 * dt * uold[M] / dx
            unew = solve(A, uold.transpose()).transpose()
            uold = unew
            u[0:M + 1, n] = unew

    return u, x, t