from collections import namedtuple
import os
import sys
import unittest
from xml.etree import ElementTree


OUTPUT_TAG = 'output'
ATOMIC_STRUCTURE_TAG = OUTPUT_TAG + '/atomic_structure'
ATOMIC_POSITIONS_TAG = ATOMIC_STRUCTURE_TAG + '/atomic_positions'
CELL_VECTORS_TAG = ATOMIC_STRUCTURE_TAG + '/cell/a%d'

TOT_ENERGY_TAG = OUTPUT_TAG + '/total_energy/etot'
ECUTWFC_TAG = OUTPUT_TAG + '/basis_set/ecutwfc'
ECUTRHO_TAG = OUTPUT_TAG + '/basis_set/ecutrho'
DFT_FUNCT_TAG = OUTPUT_TAG + '/dft/functional'
BAND_STRUCTURE_TAG = OUTPUT_TAG + '/band_structure'
NBND_TAG = BAND_STRUCTURE_TAG + '/nbnd'
NKS_TAG = BAND_STRUCTURE_TAG + '/nks'
NBND_UP_TAG = BAND_STRUCTURE_TAG + '/nbnd_up'
NBND_DW_TAG = BAND_STRUCTURE_TAG + '/nbnd_dw'
EFERMI_TAG = BAND_STRUCTURE_TAG + '/fermi_energy'
HOMO_TAG = BAND_STRUCTURE_TAG + '/highestOccupiedLevel'
KS_ENERGIES_TAG = BAND_STRUCTURE_TAG + '/ks_energies'

TOLS = dict()
TOLS[TOT_ENERGY_TAG] = 5
TOLS[EFERMI_TAG] = 5
TOLS[NBND_TAG] = 0
TOLS[NKS_TAG] = 0
TOLS[DFT_FUNCT_TAG] = -1 # str
TOLS[ECUTWFC_TAG] = 2
TOLS[ECUTRHO_TAG] = 2

XMLValidation = namedtuple('XMLValidation', ('tag', 'attrib', 'ref_value'))

# Do not use defaultdict here, since it won't throw KeyError exceptions
VALIDATION_DICT = dict()
VALIDATION_DICT['WaterP1_0_scf_0'] = list()
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(TOT_ENERGY_TAG, False, -16.91816749))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(EFERMI_TAG, False, 0.18547017))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(NBND_TAG, False, 14))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(NKS_TAG, False, 8))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(DFT_FUNCT_TAG, False, 'PBE'))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(ECUTWFC_TAG, False, 20))
VALIDATION_DICT['WaterP1_0_scf_0'].append(XMLValidation(ECUTRHO_TAG, False, 100))


class XMLTester(unittest.TestCase):
    def _assertXML(self, xml_fn, xml_validation):
        tree = ElementTree.parse(xml_fn)
        root = tree.getroot()

        ref = xml_validation.ref_value
        tol = TOLS[xml_validation.tag]

        if tol > 0:
            type_ = float
            assert_funct = lambda val: self.assertAlmostEqual(val, ref, tol)
        elif TOLS[xml_validation.tag] == 0:
            type_ = int
            assert_funct = lambda val: self.assertEqual(val, ref)
        else:
            type_ = str
            assert_funct = lambda val: self.assertEqual(val, ref)

        if xml_validation.attrib:
            val = type_(root.find(xml_validation.tag).attrib[xml_validation.attrib].strip())
        else:
            val = type_(root.find(xml_validation.tag).text.strip())

        assert_funct(val)

    def testFiles(self):
        inp_fn, ext = os.path.basename(sys.argv[1]).split('.')
        xml_fn = os.path.join(inp_fn + '.save', 'data-file-schema.xml')
        for xml_validator in VALIDATION_DICT[inp_fn]:
            self._assertXML(xml_fn, xml_validator)


if __name__ == '__main__':
    if sys.version_info < (2, 7, 0):
        sys.stderr.write("You need python 2.7 or later to run this program\n")
        sys.exit(1)
    # This is needed to run tests from the command line:
    # https://tellthemuserstories.wordpress.com/2013/01/19/python-unittest-command-line-argument/amp/
    suite = unittest.TestLoader().loadTestsFromTestCase(XMLTester)
    unittest.TextTestRunner().run(suite)
