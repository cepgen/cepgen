class PrintHelper(object):
    '''Helper class for the pretty-printing of configuration parameters'''
    _indent = 0
    _indent_size = 4
    def indent(self):
        '''Move to the next indentation block'''
        self._indent += self._indent_size
    def unindent(self):
        '''Go up to the previous indentation block'''
        self._indent -= self._indent_size
    def indentation(self):
        '''Current indentation level'''
        return ' '*self._indent

class Parameters(dict):
    '''A raw list of steering parameters'''
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def __init__(self, *args, **kwargs):
        """Construct a parameters set from dictionary arguments

        Args:
            args: list of arguments

        Kwargs:
            kwargs: list of named (keyworded) arguments

        Examples:
            >>> print dict(Parameters(foo = 'bar', foo2 = 42))
            {'foo': 'bar', 'foo2': 42}
        """
        self.update(*args, **kwargs)
        super(Parameters, self).__init__(*args, **kwargs)
    def __deepcopy__(self, memo):
        '''Override the default dict deep copy operator'''
        from copy import deepcopy
        return Parameters([(deepcopy(k, memo), deepcopy(v, memo)) for k, v in self.items()])
    def dump(self, printer=PrintHelper()):
        '''Human-readable dump of this object'''
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
        '''Human-readable version of this object'''
        return self.dump()
    def clone(self, *args, **kwargs):
        '''Return a deep copy of this object'''
        from copy import deepcopy
        out = deepcopy(self)
        for k in kwargs:
            out[k] = kwargs.get(k)
        return type(self)(out)
    def load(self, mod):
        '''Extend this object by an include'''
        mod = mod.replace('/', '.')
        module = __import__(mod)
        self.extend(sys.modules[mod])

class Module(Parameters):
    """A named parameters set to steer a generic module

    Attributes:
        mod_name: Name of this module
    """

    def __init__(self, name, *args, **kwargs):
        """Construct a module parameters set from dictionary arguments

        Args:
            name: module name
            args: list of arguments

        Kwargs:
            kwargs: list of named (keyworded) arguments

        Examples:
            >>> print dict(Module('module1', foo = 'bar', foo2 = 42))
            {'foo': 'bar', 'foo2': 42, 'mod_name': 'module1'}
        """
        super(Module, self).__init__(*args, **kwargs)
        self.mod_name = name
    def __len__(self):
        '''Number of keys handled'''
        return dict.__len__(self)-1 # discard the name key
    def dump(self, printer=PrintHelper()):
        '''Human-readable dump of this object'''
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
        '''Human-readable version of this object'''
        return self.dump()
    def clone(self, name, **kwargs):
        '''Return a deep copy of this object'''
        out = Parameters(self).clone(**kwargs)
        out.mod_name = name
        return out

class Sequence(Parameters):
    MODULE = object()
    def __init__(self, *args):
        params = Parameters()
        for mod in args:
            params[mod.mod_name] = mod
        super(Sequence, self).__init__(params)
    def __delitem__(self, index):
        self[index] = self.MODULE
    def __iter__(self):
        return (item for item in super().__iter__() if item is not self.MODULE)
    def __eq__(self, other):
        if isinstance(other, Sequence):
            return all(x == y for x, y in zip(self, other))
        return super().__eq__(other)

if __name__ == '__main__':
    import unittest
    class TestTypes(unittest.TestCase):
        '''A small collection of tests for our new types'''
        def testModules(self):
            '''Test the Module object'''
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
            '''Test the Parameters object'''
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

