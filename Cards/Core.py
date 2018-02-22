#!/usr/bin/env python

class Parameters(dict):
    '''A raw list of steering parameters'''
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

class Module(Parameters):
    '''A named parameters set to steer a generic module'''
    def __init__(self, m_name, *args, **kwargs):
        super(Module, self).__init__(*args, **kwargs)
        self.mod_name = m_name

class StructureFunctions:
    '''Types of structure functions supported'''
    Electron            = 1
    ElasticProton       = 2
    SuriYennie          = 11
    SzczurekUleshchenko = 12
    BlockDurandHa       = 13
    FioreBrasse         = 101
    ChristyBosted       = 102
    CLAS                = 103
    ALLM91              = 201
    ALLM97              = 202
    GD07p               = 203
    GD11p               = 204
    MSTWgrid            = 205
    LUXlike             = 301

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

    unittest.main()
