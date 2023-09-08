import vegas

def integrate(f, num_dim: int, num_iter: int, num_warmup: int, num_calls: int):
    integ = vegas.Integrator(num_dim * [(0., 1.)])
    def f_pyarr(vars):
        return f(vars.tolist())
    integ(f_pyarr, nitn=num_iter, neval=num_warmup)
    res = integ(f_pyarr, nitn=num_iter, neval=num_calls)
    return (res.mean, res.sdev)

if __name__ == '__main__':
    def f(x):
        return x[0]**2 + x[1]**2
    print(integrate(f, 2, 10, 1000, 1000))
