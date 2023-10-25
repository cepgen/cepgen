import torch
import numpy as np
from torchquad import MonteCarlo, set_up_backend

set_up_backend('torch', data_type='float32')

def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int, limits: list[tuple[float]]=[]):
    limits = limits if len(limits) > 0 else num_dim * [(0., 1.)]

    def func(xarr):
        return torch.from_numpy(np.array([f([float(x) for x in xvals]) for xvals in xarr.numpy()]))

    mc = MonteCarlo()
    res = mc.integrate(func, dim=num_dim, N=num_calls, integration_domain=limits, backend='torch')
    return (float(res.numpy()), 1.)

if __name__ == '__main__':
    import math
    print(integrate(lambda x: x[0]**2 + x[1]**2, 2, 10, 1000, 1000))
    print(integrate(lambda x: math.sin(x[0]), 1, 10, 1000, 1000, [(0, math.pi)]))
