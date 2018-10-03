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

