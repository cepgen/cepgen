#!/usr/bin/env python

class PrintHelper(object):
    _indent = 0
    _indent_size = 4
    def indent(self):
        self._indent += self._indent_size
    def unindent(self):
        self._indent -= self._indent_size
    def indentation(self):
        return ' '*self._indent

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
    def dump(self, printer=PrintHelper()):
        out = self.__class__.__name__+'(\n'
        for k, v in self.items():
            printer.indent()
            out += ('%s%s = ' % (printer.indentation(), k))
            if v.__class__.__name__ not in ['Parameters', 'Module']:
                out += v.__repr__()
            else:
                out += v.dump(printer)
            out += ',\n'
            printer.unindent()
        out += printer.indentation()+')'
        return out
    def __repr__(self):
        return self.dump()
    def clone(self, *args, **kwargs):
        from copy import deepcopy
        out = deepcopy(self)
        for k in kwargs:
            out[k] = kwargs.get(k)
        return type(self)(out)
    def load(self, mod):
        mod = mod.replace('/', '.')
        module = __import__(mod)
        self.extend(sys.modules[mod])

class Module(Parameters):
    '''A named parameters set to steer a generic module'''
    def __init__(self, name, *args, **kwargs):
        super(Module, self).__init__(*args, **kwargs)
        self.mod_name = name
    def __len__(self):
        return dict.__len__(self)-1 # discard the name key
    def dump(self, printer=PrintHelper()):
        out = self.__class__.__name__+'('+self.mod_name.__repr__()+'\n'
        mod_repr = self.clone('')
        mod_repr.pop('mod_name', None)
        for k, v in mod_repr.items():
            printer.indent()
            out += ('%s%s = ' % (printer.indentation(), k))
            if v.__class__.__name__ not in ['Parameters', 'Module']:
                out += v.__repr__()
            else:
                out += v.dump(printer)
            out += ',\n'
            printer.unindent()
        out += printer.indentation()+')'
        return out
    def __repr__(self):
        return self.dump()
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
    class PDFMode:
        AllQuarks     = 0
        ValenceQuarks = 1
        SeaQuarks     = 2
    '''Types of structure functions supported'''
    Electron            = Parameters(id=1)
    ElasticProton       = Parameters(id=2)
    SuriYennie          = Parameters(id=11)
    SzczurekUleshchenko = Parameters(id=12)
    BlockDurandHa       = Parameters(id=13)
    FioreBrasse         = Parameters(id=101)
    ChristyBosted       = Parameters(id=102)
    CLAS                = Parameters(id=103)
    ALLM91              = Parameters(id=201)
    ALLM97              = Parameters(id=202)
    GD07p               = Parameters(id=203)
    GD11p               = Parameters(id=204)
    MSTWgrid = Parameters(
        id = 205,
        gridPath = 'External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat',
    )
    LUXlike = Parameters(
        id = 301,
        #Q2cut = 10.,
        #W2limits = (4.,1.),
        continuumSF = GD11p,
        resonancesSF = ChristyBosted,
    )
    LHAPDF = Parameters(
        id = 401,
        pdfSet = 'LUXqed17_plus_PDF4LHC15_nnlo_100',
        numFlavours = 4,
        mode = PDFMode.AllQuarks,
    )

class KTFlux:
    PhotonElastic         = 0
    PhotonInelastic       = 1
    PhotonInelasticBudnev = 11

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

