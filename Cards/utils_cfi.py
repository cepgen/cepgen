class StructureFunctions(tuple):
    def __new__(self, name, variant=''):
        return tuple.__new__(StructureFunctions, (name, variant))

class Parameters(dict):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

class ProcessMode:
    ElectronProton = 0
    ElasticElastic = 1
    ElasticInelastic = 2
    InelasticElastic = 3
    InelasticInelastic = 4
    ProtonElectron = 5
    ElectronElectron = 6
