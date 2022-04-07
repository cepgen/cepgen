from .containers_cfi import Parameters

vermaserenCuts = Parameters(
    pt = (0.5,),
    energy = (1.,),
    eta = (-3.131, 3.131), # 5 < theta < 175 deg
    mx = (1.07, 320.),
    q2 = (0., 10000.)
)
