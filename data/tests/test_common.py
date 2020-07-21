from unittest import TestCase
from data.common import *

class SequenceSiteToVectorTests(TestCase):

    def test_can_get_features_from_sequence(self):
        sequence = "xxxxXxxxxxxxxXxx"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8
        })
        sequence = "xxxxXxxxxxxxxXxxXXxxxXxxxxxxX"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8, "gap2": 2, "gap3": 0, "gap4": 3, "gap5": 6
        })



class FamilySplittingTests(TestCase):

    def test_can_split_single_residue_family(self):
        self.assertEqual(split_family("H3"), [["H", 3]])


    def test_can_split_multi_residue_family(self):
        self.assertEqual(split_family("H3C10"), [["H", 3], ["C", 10]])