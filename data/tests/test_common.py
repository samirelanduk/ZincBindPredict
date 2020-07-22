from atomium import Residue, Atom, Chain
from unittest import TestCase
from unittest.mock import patch
from data.common import *

class SequenceSiteToVectorTests(TestCase):

    @patch("data.common.average_hydrophobicity")
    def test_can_get_features_from_sequence(self, mock_hydro):
        mock_hydro.side_effect = [1.1, 0.5]
        sequence = "xxxxXxxxxxxxxXxx"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8, "hydrophobicity": 1.1
        })
        mock_hydro.assert_called_with(sequence)
        sequence = "xxxxXxxxxxxxxXxxXXxxxXxxxxxxX"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8, "gap2": 2, "gap3": 0, "gap4": 3, "gap5": 6,
            "hydrophobicity": 0.5
        })
        mock_hydro.assert_called_with(sequence)



class FamilySplittingTests(TestCase):

    def test_can_split_single_residue_family(self):
        self.assertEqual(split_family("H3"), [["H", 3]])


    def test_can_split_multi_residue_family(self):
        self.assertEqual(split_family("H3C10"), [["H", 3], ["C", 10]])



class AverageHydrophobicityTests(TestCase):

    def test_average_hydrophobicity_one_residue(self):
        self.assertEqual(average_hydrophobicity("xaCwx"), -0.84)
    

    def test_average_hydrophobicity_multiple_residues(self):
        self.assertEqual(average_hydrophobicity("xaCwxdHt"), -0.078)
    

    def test_average_hydrophobicity_when_residues_at_end(self):
        self.assertEqual(average_hydrophobicity("xaCwxdH"), -0.15)
        self.assertEqual(average_hydrophobicity("CwxdHt"), -0.16)
        self.assertEqual(average_hydrophobicity("CwxdH"), -0.31)
    

    def test_average_hydrophobicity_unknown_residue(self):
        self.assertEqual(average_hydrophobicity("xzCwx"), -1.85)



class StructureSiteToVectorTests(TestCase):

    def setUp(self):
        self.res1 = Residue(
         Atom("C", 0, -2, 0, 1, "CA", 0, 0, []), Atom("C", 0, -1, 0, 1, "CB", 0, 0, [])
        )
        self.res2 = Residue(
         Atom("C", 0, 2, 0, 1, "CA", 0, 0, []), Atom("C", 0, 1, 0, 1, "CB", 0, 0, [])
        )
        self.res3 = Residue(
         Atom("C", -2, 0, 0, 1, "CA", 0, 0, []), Atom("C", -1, 0, 0, 1, "CB", 0, 0, [])
        )
        self.res4 = Residue(
         Atom("C", 2, 0, 0, 1, "CA", 0, 0, []), Atom("C", 1, 0, 0, 1, "CB", 0, 0, [])
        )
        Chain(
         self.res1, self.res2, self.res3, self.res4,
         helices=[[self.res2, self.res3]], strands=[[self.res4]]
        )
        

    def test_can_get_vector_dict(self):
        sample = structure_family_site_to_vector((self.res1, self.res2, self.res3, self.res4))
        self.assertEqual(sample.keys(), {
         "ca_mean", "ca_std", "ca_max", "ca_min", "cb_mean", "cb_std",
         "cb_max", "cb_min", "helix", "strand"
        })
        self.assertAlmostEqual(sample["ca_mean"], 3.218, delta=0.005)
        self.assertAlmostEqual(sample["ca_std"], 0.552, delta=0.005)
        self.assertEqual(sample["ca_max"], 4)
        self.assertAlmostEqual(sample["ca_min"], 2.828, delta=0.005)
        self.assertAlmostEqual(sample["cb_mean"], 1.609, delta=0.005)
        self.assertAlmostEqual(sample["cb_std"], 0.276, delta=0.005)
        self.assertEqual(sample["cb_max"], 2)
        self.assertAlmostEqual(sample["cb_min"], 1.414, delta=0.005)
        self.assertEqual(sample["helix"], 2)
        self.assertEqual(sample["strand"], 1)