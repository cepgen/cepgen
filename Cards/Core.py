#!/usr/bin/env python

class StructureFunctions(tuple):
    '''A set of structure functions with a small granularity'''
    def __new__(self, name, variant=''):
        return tuple.__new__(StructureFunctions, (name, variant))

class Parameters(dict):
    '''A raw list of steering parameters'''
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

class Module(dict):
    '''A named parameters set to steer a generic module'''
    def __init__(self, mname, *args, **kwargs):
        self.update(*args, **kwargs)
        self['mod_name'] = mname

class ProcessMode:
    '''Types of processes supported'''
    ElectronProton = 0
    ElasticElastic = 1
    ElasticInelastic = 2
    InelasticElastic = 3
    InelasticInelastic = 4
    ProtonElectron = 5
    ElectronElectron = 6
