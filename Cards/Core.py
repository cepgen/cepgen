#!/usr/bin/env python

class Parameters(dict):
    '''A raw list of steering parameters'''
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

class Module(Parameters):
    '''A named parameters set to steer a generic module'''
    def __init__(self, mname, *args, **kwargs):
        super(Module, self).__init__(*args, **kwargs)
        self.mod_name = mname

class StructureFunctions(tuple):
    '''A set of structure functions with a small granularity'''
    def __new__(self, name, variant=''):
        return tuple.__new__(StructureFunctions, (name, variant))
    def name(self):
        return self.__getitem__(0)
    def variant(self):
        return self.__getitem__(1)

class ProcessMode:
    '''Types of processes supported'''
    ElectronProton = 0
    ElasticElastic = 1
    ElasticInelastic = 2
    InelasticElastic = 3
    InelasticInelastic = 4
    ProtonElectron = 5
    ElectronElectron = 6

if __name__ == '__main__':
    import unittest
    class TestTypes(unittest.TestCase):
        def testModules(self):
            mod = Module('empty')
            self.assertEqual(len(mod), 1)
            mod.param1 = 'foo'
            self.assertEqual(len(mod), 2)
        def testParameters(self):
            params = Parameters(
                first = 'foo',
                second = 'bar',
                third = 42,
                fourth = (1, 2),
            )
            self.assertEqual(len(params), 4)
            self.assertEqual(params.first, params['first'])
            self.assertEqual(params['second'], 'bar')
            self.assertTrue(int(params.third) == params.third)
            self.assertEqual(len(params.fourth), 2)
        def testStructureFunctions(self):
            sf = StructureFunctions('my_sf', 'nope')
            self.assertEqual(sf.name(), 'my_sf')
            self.assertEqual(sf.variant(), 'nope')

    unittest.main()
