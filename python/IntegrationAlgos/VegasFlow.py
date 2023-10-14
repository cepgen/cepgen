from vegasflow import VegasFlow, run_eager
run_eager()
import tensorflow as tf
import numpy as np


def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int, limits: list[tuple[float]]=[]):
    limits = limits if len(limits) > 0 else num_dim * [(0., 1.)]
    limits_len = [lim[1] - lim[0] for lim in limits]
    jacob = np.prod(limits_len)
    print(limits,limits_len)
    def xnorm(xarr):
        return [limits[i][0] + limits_len[i] * xarr[i] for i in range(len(xarr))]

    vegas = VegasFlow(num_dim, num_calls, main_dimension=1)

    @tf.py_function(Tout=tf.float64)
    def func(xarr):
        return tf.constant(np.array([f(xnorm(x)) for x in xarr.numpy()]))

    vegas.compile(func)
    res, unc = vegas.run_integration(2)
    return (res * jacob, unc * jacob)

if __name__ == '__main__':
    import math
    print(integrate(lambda x: x[0]**2 + x[1]**2, 2, 10, 1000, 1000))
    print(integrate(lambda x: math.sin(x[0]), 1, 10, 1000, 1000, [(0, math.pi)]))
