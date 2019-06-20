from unittest import TestCase
from unittest.mock import patch, Mock
from core.utilities import *

class FamilySplittingTests(TestCase):

    def test_can_split_simple_family(self):
        self.assertEqual(split_family("H3"), ["H3"])
        self.assertEqual(split_family("P13"), ["P13"])
    

    def test_can_split_simple_families_with_different_residues(self):
        self.assertEqual(split_family("C4H2"), ["C4", "H2"])
        self.assertEqual(split_family("A12B6C176D1"), ["A12", "B6", "C176", "D1"])



class ModelToResidueCombinationsTests(TestCase):

    def setUp(self):
        self.patch1 = patch("core.utilities.split_family")
        self.mock_split = self.patch1.start()
        self.mock_split.side_effect = lambda family: [family]
        self.model = Mock()
        self.model.residues.side_effect = lambda code: []
    

    def tearDown(self):
        self.patch1.stop()


    def test_can_handle_no_matching_residues(self):
        self.model.residues.side_effect = lambda code: []
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, ())
    

    def test_can_handle_insufficient_matching_residues(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2"]
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, ())
    

    def test_can_return_single_combination(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2", "R3"]
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, (("R1", "R2", "R3"),))
    

    def test_can_return_many_combinations(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2", "R3", "R4", "R5"]
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, (
         ("R1", "R2", "R3"), ("R1", "R2", "R4"), ("R1", "R2", "R5"),
         ("R1", "R3", "R4"), ("R1", "R3", "R5"), ("R1", "R4", "R5"),
         ("R2", "R3", "R4"), ("R2", "R3", "R5"), ("R2", "R4", "R5"),
         ("R3", "R4", "R5")
        ))
        combos = model_to_residue_combos(self.model, "C4")
        self.mock_split.assert_called_with("C4")
        self.model.residues.assert_called_with(code="C")
        self.assertEqual(combos, (
         ("R1", "R2", "R3", "R4"), ("R1", "R2", "R3", "R5"), ("R1", "R2", "R4", "R5"),
         ("R1", "R3", "R4", "R5"), ("R2", "R3", "R4", "R5")
        ))
    

    def test_can_return_combinations_from_different_subfamilies(self):
        self.mock_split.side_effect = lambda family: ["H2", "C3"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3", "C4"]]
        combos = model_to_residue_combos(self.model, "H2C3")
        self.mock_split.assert_called_with("H2C3")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.assertEqual(combos, (
         ("H1", "H2", "C1", "C2", "C3"), ("H1", "H2", "C1", "C2", "C4"),
         ("H1", "H2", "C1", "C3", "C4"), ("H1", "H2", "C2", "C3", "C4"),
         ("H1", "H3", "C1", "C2", "C3"), ("H1", "H3", "C1", "C2", "C4"),
         ("H1", "H3", "C1", "C3", "C4"), ("H1", "H3", "C2", "C3", "C4"),
         ("H2", "H3", "C1", "C2", "C3"), ("H2", "H3", "C1", "C2", "C4"),
         ("H2", "H3", "C1", "C3", "C4"), ("H2", "H3", "C2", "C3", "C4")
        ))

        self.mock_split.side_effect = lambda family: ["H3", "C1"]
        self.model.residues.side_effect = [["H1", "H2", "H3", "H4"], ["C1", "C2"]]
        combos = model_to_residue_combos(self.model, "H3C1")
        self.mock_split.assert_called_with("H3C1")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.assertEqual(combos, (
         ("H1", "H2", "H3", "C1"), ("H1", "H2", "H3", "C2"),
         ("H1", "H2", "H4", "C1"), ("H1", "H2", "H4", "C2"),
         ("H1", "H3", "H4", "C1"), ("H1", "H3", "H4", "C2"),
         ("H2", "H3", "H4", "C1"), ("H2", "H3", "H4", "C2")
        ))

        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], ["E1", "E2", "E3"]]
        combos = model_to_residue_combos(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, (
         ("H1", "H2", "C1", "C2", "E1", "E2"), ("H1", "H2", "C1", "C2", "E1", "E3"),
         ("H1", "H2", "C1", "C2", "E2", "E3"), ("H1", "H2", "C1", "C3", "E1", "E2"),
         ("H1", "H2", "C1", "C3", "E1", "E3"), ("H1", "H2", "C1", "C3", "E2", "E3"),
         ("H1", "H2", "C2", "C3", "E1", "E2"), ("H1", "H2", "C2", "C3", "E1", "E3"),
         ("H1", "H2", "C2", "C3", "E2", "E3"), ("H1", "H3", "C1", "C2", "E1", "E2"),
         ("H1", "H3", "C1", "C2", "E1", "E3"), ("H1", "H3", "C1", "C2", "E2", "E3"),
         ("H1", "H3", "C1", "C3", "E1", "E2"), ("H1", "H3", "C1", "C3", "E1", "E3"),
         ("H1", "H3", "C1", "C3", "E2", "E3"), ("H1", "H3", "C2", "C3", "E1", "E2"),
         ("H1", "H3", "C2", "C3", "E1", "E3"), ("H1", "H3", "C2", "C3", "E2", "E3"),
         ("H2", "H3", "C1", "C2", "E1", "E2"), ("H2", "H3", "C1", "C2", "E1", "E3"),
         ("H2", "H3", "C1", "C2", "E2", "E3"), ("H2", "H3", "C1", "C3", "E1", "E2"),
         ("H2", "H3", "C1", "C3", "E1", "E3"), ("H2", "H3", "C1", "C3", "E2", "E3"),
         ("H2", "H3", "C2", "C3", "E1", "E2"), ("H2", "H3", "C2", "C3", "E1", "E3"),
         ("H2", "H3", "C2", "C3", "E2", "E3")
        ))
    

    def test_can_handle_insufficient_residues_in_one_subfamily(self):
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], []]
        combos = model_to_residue_combos(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, ())



class ResiduesToSampleTests(TestCase):

    def setUp(self):
        for r in range(1, 5):
            res = Mock()
            setattr(self, f"res{r}", res)
            ca = Mock()
            res.atom.return_value = ca
            ca.distance_to.side_effect = [1, 1.5, 2]
        

    def test_can_get_sample_dict(self):
        sample = residues_to_sample((self.res1, self.res2, self.res3, self.res4))
        self.res1.atom.assert_called_with(name="CA")
        self.res2.atom.assert_called_with(name="CA")
        self.res3.atom.assert_called_with(name="CA")
        self.res4.atom.assert_called_with(name="CA")
        self.assertEqual(sample.keys(), {"mean_ca", "ca_std"})
        self.assertEqual(sample["mean_ca"], 4 / 3)
        self.assertAlmostEqual(sample["ca_std"], 0.372, delta=0.005)
        