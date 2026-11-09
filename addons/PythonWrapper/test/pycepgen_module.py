#!/usr/bin/env python3
"""Tests for the pycepgen Python bindings"""

import math
import pycepgen as cg
import unittest


class TestPDG(unittest.TestCase):
    def test_mass(self):
        self.assertTrue(math.isclose(cg.PDG.mass(13), 0.1056583755, rel_tol=1.e-3))
    def test_charge(self):
        self.assertTrue(math.isclose(cg.PDG.charge(11), -1., rel_tol=1.e-6))
    def test_colours(self):
        self.assertEqual(cg.PDG.colours(13), 1.)


class TestParametersConversion(unittest.TestCase):
    def test_structure_functions_attributes(self):
        plain = cg.StructureFunctionsFactory.build('SuriYennie')
        reference = plain.F2(0.1, 10.)
        self.assertTrue(isinstance(plain.parameters, dict))
        self.assertEqual(plain.parameters['mod_name'], 'SuriYennie')
        # build from a parameters dictionary, with or without the explicit module name
        from_dict = cg.StructureFunctionsFactory.build({'mod_name': 'SuriYennie'})
        from_name_dict = cg.StructureFunctionsFactory.build('SuriYennie', {})
        self.assertEqual(from_dict.F2(0.1, 10.), reference)
        self.assertEqual(from_name_dict.F2(0.1, 10.), reference)
    def test_plain_parameters(self):
        # re-parse parameters
        plain = cg.StructureFunctionsFactory.build('SuriYennie')
        plain.parameters = {'mod_name': 'SuriYennie', 'testInt': 1, 'testRange': (0.5, 1.5), 'testList': [{'a': 1}]}
        self.assertEqual(plain.parameters['testInt'], 1)
        self.assertEqual(plain.parameters['testRange'], (0.5, 1.5))
    def test_invalid_arguments(self):
        with self.assertRaises(TypeError):
            cg.StructureFunctionsFactory.build([])


class TestStructureFunctions(unittest.TestCase):
    _sf_by_name = cg.StructureFunctionsFactory.build('SuriYennie')
    _sf_by_index = cg.StructureFunctionsFactory.build(11)
    def test_modules_content(self):
        self.assertTrue('SuriYennie' in cg.StructureFunctionsFactory.modules)
    def test_sfs_by_name(self):
        self._test_sfs(self._sf_by_name)
    def test_sfs_by_index(self):
        self._test_sfs(self._sf_by_index)
    def test_modules_usage(self):
        self.assertEqual(self._sf_by_name.F2(0.1, 10.), self._sf_by_index.F2(0.1, 10.))
    def _test_sfs(self, sfs):
        self.assertTrue(isinstance(sfs, cg._StructureFunctions))
        self.assertTrue(isinstance(sfs, cg._Steerable))
        self.assertEqual(sfs.name, 'SuriYennie')
        self.assertTrue(sfs.F2(0.1, 10.) > 0.)
        self.assertTrue(math.isfinite(sfs.FL(0.1, 10.)))
        self.assertTrue(math.isfinite(sfs.F1(0.1, 10.)))


class TestSigmaRatios(unittest.TestCase):
    def test_sigma_ratios(self):
        index = 1
        ratio, uncertainty = cg.SigmaRatiosFactory.build(index)(0.1, 10.)
        self.assertTrue(isinstance(cg.SigmaRatiosFactory.build(index), cg._SigmaRatio))
        self.assertTrue(math.isfinite(ratio))
        self.assertTrue(math.isfinite(uncertainty))
        self.assertTrue(cg.SigmaRatiosFactory.modules)


class TestFormFactors(unittest.TestCase):
    def test_has_form_factors(self):
        self.assertTrue(len(cg.FormFactorsFactory.modules) > 0)
    def test_form_factors(self):
        form_factors = cg.FormFactorsFactory.build('StandardDipole')
        values = form_factors(1.)
        self.assertTrue(isinstance(values, cg.FormFactors))
        self.assertTrue(math.isfinite(values.FE))
        self.assertTrue(values.GE > 0.)
        self.assertTrue(math.isfinite(values.FM))
        self.assertTrue(values.GM > 0.)


class TestPartonFluxes(unittest.TestCase):
    def test_has_kt_parton_fluxes(self):
        self.assertTrue(len(cg.KTFluxFactory.modules) > 0)
    def test_kt_parton_fluxes(self):
        kt = cg.KTFluxFactory.build('BudnevInelastic')
        self.assertTrue(isinstance(kt, cg._KTFlux))
        self.assertTrue(kt.ktFactorised)
        self.assertTrue(kt() is kt)
        self.assertTrue(math.isfinite(kt.fluxQ2(0.1, 1., 2.)))
        self.assertTrue(math.isfinite(kt.fluxMX2(0.1, 1., 2.)))
    def test_has_coll_parton_fluxes(self):
        self.assertTrue(len(cg.CollinearFluxFactory.modules) > 0)
    def test_coll_parton_fluxes(self):
        collinear = cg.CollinearFluxFactory.build('EPAFlux')
        self.assertTrue(isinstance(collinear, cg._CollinearFlux))
        self.assertFalse(collinear.ktFactorised)
        self.assertTrue(collinear() is collinear)
        self.assertTrue(math.isfinite(collinear.fluxMX2(0.1, 0.)))
        self.assertTrue(isinstance(collinear.partonPdgId, int))


if __name__ == '__main__':
    unittest.main()
