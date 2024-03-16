##
# \file
# \ingroup python integration
#
# MCint integration algorithm interface

import mcint
import random
import numpy as np


def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int, limits: list[tuple[float]]=[]):
    limits = limits if len(limits) > 0 else num_dim * [(0., 1.)]
    jacob = np.prod([lim[1] - lim[0] for lim in limits])
    def sampler():
        while True:
            yield [random.uniform(lim[0], lim[1]) for lim in limits]
    return mcint.integrate(f, sampler(), measure=jacob, n=num_calls)


if __name__ == '__main__':
    import math
    print(integrate(lambda x: x[0]**2 + x[1]**2, 2, 10, 1000, 1000))
    print(integrate(lambda x: math.sin(x[0]), 1, 10, 1000, 1000, [(0, math.pi)]))
