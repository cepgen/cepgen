#!/usr/bin/env python

'''@package cepgen
A collection of tools for Python steering cards definition
'''

from containers_cfi import Module, Parameters
#--- physics-level includes
from StructureFunctions_cfi import StructureFunctions
from ProcessMode_cfi import ProcessMode

class Logging:
    '''Logging verbosity'''
    Nothing         = 0
    Error           = 1
    Warning         = 2
    Information     = 3
    Debug           = 4
    DebugInsideLoop = 5

if __name__ == '__main__':
    import unittest
    class TestTypes(unittest.TestCase):
        def testModules(self):
            mod = Module('empty')
            self.assertEqual(len(mod), 0)
            mod.param1 = 'foo'
            self.assertEqual(len(mod), 1)
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
            self.assertEqual(params.third, 43)
            self.assertEqual(params_copy.third, 42)

    unittest.main()

