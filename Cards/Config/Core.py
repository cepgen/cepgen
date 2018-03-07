#!/usr/bin/env python

class Parameters(dict):
    '''A raw list of steering parameters'''
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        super(Parameters, self).__init__(*args, **kwargs)
    def __deepcopy__(self, memo):
        from copy import deepcopy
        return Parameters([(deepcopy(k, memo), deepcopy(v, memo)) for k, v in self.items()])
    def clone(self, *args, **kwargs):
        from copy import deepcopy
        out = deepcopy(self)
        for k in kwargs:
            out[k] = kwargs.get(k)
        return type(self)(out)


class Module(Parameters):
    '''A named parameters set to steer a generic module'''
    def __init__(self, name, *args, **kwargs):
        super(Module, self).__init__(*args, **kwargs)
        self.mod_name = name
    def clone(self, name, **kwargs):
        out = Parameters(self).clone(**kwargs)
        out.mod_name = name
        return out

class Logging:
    '''Logging verbosity'''
    Nothing         = 0
    Error           = 1
    Warning         = 2
    Information     = 3
    Debug           = 4
    DebugInsideLoop = 5

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
            # playing with modules clones
            mod_copy = mod.clone('notEmpty', param1 = 'boo', param2 = 'bar')
            self.assertEqual(mod.param1, 'foo')
            self.assertEqual(mod_copy.param1, 'boo')
            self.assertEqual(mod_copy.param2, 'bar')
            self.assertEqual(mod.param1+mod_copy.param2, 'foobar')
        def testParameters(self):
            params = Parameters(
                first = 'foo',
                second = 'bar',
                third = 42,
                fourth = (1, 2),
            )
            params_copy = params.clone(
                second = 'bak',
            )
            self.assertEqual(len(params), 4)
            self.assertEqual(params.first, params['first'])
            self.assertEqual(params['second'], 'bar')
            self.assertTrue(int(params.third) == params.third)
            self.assertEqual(len(params.fourth), 2)
            self.assertEqual(params.second, 'bar')
            # playing with parameters clones
            self.assertEqual(params_copy.second, 'bak')
            # check that the clone does not change value if the origin does
            # (i.e. we indeed have a deep copy and not a shallow one...)
            params.third = 43
            self.assertEqual(params.third, 43 )
            self.assertEqual(params_copy.third, 42 )

    unittest.main()
