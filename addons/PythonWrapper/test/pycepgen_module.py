#!/usr/bin/env python3
"""Tests for the pycepgen Python bindings"""

import math
import pycepgen as cg


def test_pdg():
    assert math.isclose(cg.PDG.mass(13), 0.1056583755, rel_tol=1.e-3)
    assert math.isclose(cg.PDG.charge(11), -1., rel_tol=1.e-6) or math.isclose(cg.PDG.charge(11), 1., rel_tol=1.e-6)
    assert cg.PDG.colours(13) == 1.


def test_structure_functions():
    assert 'SuriYennie' in cg.StructureFunctionsFactory.modules
    by_name, by_index = cg.StructureFunctionsFactory.build('SuriYennie'), cg.StructureFunctionsFactory.build(11)
    assert isinstance(by_name, cg._StructureFunctions) and isinstance(by_name, cg._Steerable)
    assert by_name.name == 'SuriYennie' and by_index.name == 'SuriYennie'
    assert by_name.F2(0.1, 10.) == by_index.F2(0.1, 10.) > 0.
    assert math.isfinite(by_name.FL(0.1, 10.)) and math.isfinite(by_name.F1(0.1, 10.))
    assert 'SuriYennie' in cg.StructureFunctionsFactory.describe('SuriYennie')


def test_parameters_conversion():
    plain = cg.StructureFunctionsFactory.build('SuriYennie')
    reference = plain.F2(0.1, 10.)
    assert isinstance(plain.parameters, dict) and plain.parameters['mod_name'] == 'SuriYennie'
    # build from a parameters dictionary, with or without the explicit module name
    from_dict = cg.StructureFunctionsFactory.build({'mod_name': 'SuriYennie'})
    from_name_dict = cg.StructureFunctionsFactory.build('SuriYennie', {})
    assert from_dict.F2(0.1, 10.) == reference and from_name_dict.F2(0.1, 10.) == reference
    plain.parameters = {'mod_name': 'SuriYennie', 'testInt': 1, 'testRange': (0.5, 1.5), 'testList': [{'a': 1}]}
    assert plain.parameters['testInt'] == 1
    assert plain.parameters['testRange'] == (0.5, 1.5)


def test_sigma_ratios():
    index = 1
    ratio, uncertainty = cg.SigmaRatiosFactory.build(index)(0.1, 10.)
    assert isinstance(cg.SigmaRatiosFactory.build(index), cg._SigmaRatio)
    assert math.isfinite(ratio) and math.isfinite(uncertainty)
    assert cg.SigmaRatiosFactory.modules


def test_form_factors():
    assert cg.FormFactorsFactory.modules
    form_factors = cg.FormFactorsFactory.build('StandardDipole')
    values = form_factors(1.)
    assert isinstance(values, cg.FormFactors)
    assert math.isfinite(values.FE) and math.isfinite(values.FM)
    assert values.GE > 0. and values.GM > 0.


def test_parton_fluxes():
    assert cg.KTFluxFactory.modules and cg.CollinearFluxFactory.modules
    collinear = cg.CollinearFluxFactory.build('EPAFlux')
    assert isinstance(collinear, cg._CollinearFlux) and not collinear.ktFactorised
    assert collinear() is collinear
    assert math.isfinite(collinear.fluxMX2(0.1, 0.)) and isinstance(collinear.partonPdgId, int)
    kt = cg.KTFluxFactory.build('BudnevInelastic')
    assert isinstance(kt, cg._KTFlux) and kt.ktFactorised and kt() is kt
    assert math.isfinite(kt.fluxQ2(0.1, 1., 2.)) and math.isfinite(kt.fluxMX2(0.1, 1., 2.))


def test_invalid_arguments():
    try:
        cg.StructureFunctionsFactory.build([])
    except TypeError:
        return
    raise AssertionError('a list should not have been accepted as a module definition')


if __name__ == '__main__':
    for name, test in sorted(globals().items()):
        if name.startswith('test_'):
            test()
            print(f'Passed {name}')
