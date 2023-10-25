import vegas

def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int, limits: list[tuple[float]]=[]):
    limits = limits if len(limits) > 0 else num_dim * [(0., 1.)]
    integ = vegas.Integrator(limits)
    def f_pyarr(vars):
        return f(vars.tolist())
    integ(f_pyarr, nitn=num_iter, neval=num_warmup)
    res = integ(f_pyarr, nitn=num_iter, neval=num_calls)
    return (res.mean, res.sdev)

if __name__ == '__main__':
    import math
    print(integrate(lambda x: x[0]**2 + x[1]**2, 2, 10, 1000, 1000))
    print(integrate(lambda x: math.sin(x[0]), 1, 10, 1000, 1000, [(0, math.pi)]))
