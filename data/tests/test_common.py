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