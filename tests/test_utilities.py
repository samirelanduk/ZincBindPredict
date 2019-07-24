from atomium import Residue, Atom, Chain
from unittest import TestCase
from unittest.mock import patch, Mock
from core.utilities import *

class DataFetchingTests(TestCase):

    @patch("kirjava.execute")
    def test_can_get_data(self, mock_ex):
        mock_ex.return_value = {"data": {"objects": {"edges":[
         {"node": {"id": "123", "subObjects": {"edges": [
          {"node": {"id": "1"}}, {"node": {"id": "2"}}
         ]}}},
         {"node": {"id": "1234", "subObjects": {"edges": [
          {"node": {"id": "29"}}
         ]}}}
        ]}}}
        data = fetch_data("http://", "{query {1 2}}", {1: 2})
        mock_ex.assert_called_with("http://", "{query {1 2}}", variables={1: 2})
        self.assertEqual(data, [
         {"id": "123", "subObjects": [{"id": "1"}, {"id": "2"}]},
         {"id": "1234", "subObjects": [{"id": "29"}]}
        ])


class FamilySplittingTests(TestCase):

    def test_can_split_simple_family(self):
        self.assertEqual(split_family("H3"), ["H3"])
        self.assertEqual(split_family("P13"), ["P13"])
    

    def test_can_split_simple_families_with_different_residues(self):
        self.assertEqual(split_family("C4H2"), ["C4", "H2"])
        self.assertEqual(split_family("A12B6C176D1"), ["A12", "B6", "C176", "D1"])



class ModelCombinationsCountTests(TestCase):

    def setUp(self):
        self.patch1 = patch("core.utilities.split_family")
        self.mock_split = self.patch1.start()
        self.mock_split.side_effect = lambda family: [family]
        self.model = Mock()
        self.model.residues.side_effect = lambda code: []
    

    def tearDown(self):
        self.patch1.stop()


    def test_can_count_no_matching_residues(self):
        self.model.residues.side_effect = lambda code: []
        count = count_combinations(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(count, 0)
    

    def test_can_count_insufficient_matching_residues(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2"]
        count = count_combinations(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(count, 0)
    

    def test_can_count_single_combination(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2", "R3"]
        count = count_combinations(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(count, 1)
    

    def test_can_count_many_combinations(self):
        self.model.residues.side_effect = lambda code: ["R1", "R2", "R3", "R4", "R5"]
        count = count_combinations(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(count, 10)
        count = count_combinations(self.model, "C4")
        self.mock_split.assert_called_with("C4")
        self.model.residues.assert_called_with(code="C")
        self.assertEqual(count, 5)
    

    def test_can_return_combinations_from_different_subfamilies(self):
        self.mock_split.side_effect = lambda family: ["H2", "C3"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3", "C4"]]
        count = count_combinations(self.model, "H2C3")
        self.mock_split.assert_called_with("H2C3")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.assertEqual(count, 12)

        self.mock_split.side_effect = lambda family: ["H3", "C1"]
        self.model.residues.side_effect = [["H1", "H2", "H3", "H4"], ["C1", "C2"]]
        count = count_combinations(self.model, "H3C1")
        self.mock_split.assert_called_with("H3C1")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.assertEqual(count, 8)

        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], ["E1", "E2", "E3"]]
        count = count_combinations(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(count, 27)
    

    def test_can_count_insufficient_residues_in_one_subfamily(self):
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], []]
        count = count_combinations(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(count, 0)


'''
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


    def test_can_limit_number_of_combinations_returned(self):
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], ["E1", "E2", "E3"]]
        combos = model_to_residue_combos(self.model, "H2C2E2", limit=4)
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, (
         ("H1", "H2", "C1", "C2", "E1", "E2"), ("H1", "H2", "C1", "C2", "E1", "E3"),
         ("H1", "H2", "C1", "C2", "E2", "E3"), ("H1", "H2", "C1", "C3", "E1", "E2"),
        ))
'''


class SequenceToResidueComboTests(TestCase):

    def setUp(self):
        self.patch1 = patch("core.utilities.split_family")
        self.mock_split = self.patch1.start()
        self.mock_split.side_effect = lambda family: [family]
    

    def tearDown(self):
        self.patch1.stop()


    def test_can_handle_no_matching_residues(self):
        combos = sequence_to_residue_combos("abcd", "H3")
        self.mock_split.assert_called_with("H3")
        self.assertEqual(combos, [])
    

    def test_can_handle_insufficient_matching_residues(self):
        combos = sequence_to_residue_combos("abcHdeHfgi", "H3")
        self.mock_split.assert_called_with("H3")
        self.assertEqual(combos, [])
    

    def test_can_return_single_combination(self):
        combos = sequence_to_residue_combos("abcHdeHfghi", "H3")
        self.mock_split.assert_called_with("H3")
        self.assertEqual(combos, ["abcHdeHfgHi"])
    

    def test_can_return_many_combinations(self):
        combos = sequence_to_residue_combos("abcHdeHfghihsdsd", "H3")
        self.mock_split.assert_called_with("H3")
        self.assertEqual(combos, [
         "abcHdeHfgHihsdsd", "abcHdeHfghiHsdsd", "abcHdehfgHiHsdsd", "abchdeHfgHiHsdsd"
        ])
    

    def test_can_return_combinations_from_different_subfamilies(self):
        self.mock_split.side_effect = lambda family: ["H2", "C3"]
        combos = sequence_to_residue_combos("xhxyzhxyzxhxycxyzcxyzxyzcxycxyz", "H2C3")
        self.mock_split.assert_called_with("H2C3")
        self.assertEqual(combos, [
         "xHxyzHxyzxhxyCxyzCxyzxyzCxycxyz", "xHxyzHxyzxhxyCxyzCxyzxyzcxyCxyz",
         "xHxyzHxyzxhxyCxyzcxyzxyzCxyCxyz", "xHxyzHxyzxhxycxyzCxyzxyzCxyCxyz",
         "xHxyzhxyzxHxyCxyzCxyzxyzCxycxyz", "xHxyzhxyzxHxyCxyzCxyzxyzcxyCxyz",
         "xHxyzhxyzxHxyCxyzcxyzxyzCxyCxyz", "xHxyzhxyzxHxycxyzCxyzxyzCxyCxyz",
         "xhxyzHxyzxHxyCxyzCxyzxyzCxycxyz", "xhxyzHxyzxHxyCxyzCxyzxyzcxyCxyz",
         "xhxyzHxyzxHxyCxyzcxyzxyzCxyCxyz", "xhxyzHxyzxHxycxyzCxyzxyzCxyCxyz"
        ])
    

    def test_can_handle_insufficient_residues_in_one_subfamily(self):
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        combos = sequence_to_residue_combos("xhxhxhxcxexexe", "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.assertEqual(combos, [])
    

    def test_can_limit_number_of_combinations_returned(self):
        self.mock_split.side_effect = lambda family: ["H2", "C3"]
        combos = sequence_to_residue_combos("xhxyzhxyzxhxycxyzcxyzxyzcxycxyz", "H2C3", limit=4)
        self.mock_split.assert_called_with("H2C3")
        self.assertEqual(combos, [
         "xHxyzHxyzxhxyCxyzCxyzxyzCxycxyz", "xHxyzHxyzxhxyCxyzCxyzxyzcxyCxyz",
         "xHxyzHxyzxhxyCxyzcxyzxyzCxyCxyz", "xHxyzHxyzxhxycxyzCxyzxyzCxyCxyz",
        ])



class ResiduesToSampleTests(TestCase):

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
        

    def test_can_get_sample_dict(self):
        sample = residues_to_sample((self.res1, self.res2, self.res3, self.res4), "X1")
        self.assertEqual(sample.keys(), {
         "site", "ca_mean", "ca_std", "ca_max", "ca_min", "cb_mean", "cb_std",
         "cb_max", "cb_min", "helix", "strand"
        })
        self.assertEqual(sample["site"], "X1")
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



class SequenceToSampleTests(TestCase):

    def test_can_get_sequence_dict(self):
        sample = sequence_to_sample("abcHsdHaslkdHoiuhdsHasw", "X1")
        self.assertEqual(sample.keys(), {
         "site", "min_spacer", "max_spacer", "spacer_1", "spacer_2", "spacer_3"
        })
        self.assertEqual(sample["site"], "X1")
        self.assertEqual(sample["min_spacer"], 2)
        self.assertEqual(sample["max_spacer"], 6)
        self.assertEqual(sample["spacer_1"], 2)
        self.assertEqual(sample["spacer_2"], 5)
        self.assertEqual(sample["spacer_3"], 6)
        