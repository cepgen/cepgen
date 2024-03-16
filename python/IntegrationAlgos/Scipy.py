##
# \file
# \ingroup python integration
#
# Scipy integration algorithm interface

from scipy import integrate as spint
import numpy as np


def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int, limits: list[tuple[float]]=[]):
    limits = [[lim[0], lim[1]] for lim in (limits if len(limits) > 0 else num_dim * [[0., 1.]])]
    def f_args(*args):
        return f(args)
    return spint.nquad(f_args, limits)


if __name__ == '__main__':
    import math
    print(integrate(lambda x: x[0]**2 + x[1]**2, 2, 10, 1000, 1000))
    print(integrate(lambda x: math.sin(x[0]), 1, 10, 1000, 1000, [(0, math.pi)]))
