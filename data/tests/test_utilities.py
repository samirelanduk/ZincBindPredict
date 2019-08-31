from atomium import Residue, Atom, Chain
import pandas as pd
import numpy as np
from io import StringIO
from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from data.utilities import *

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



class DataFileUpdatingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("builtins.open")
        self.mock_open = self.patch1.start()
        open_return = MagicMock()
        mock_file = Mock()
        self.mock_write = MagicMock()
        mock_file.write = self.mock_write
        open_return.__enter__.return_value = mock_file
        self.mock_writelines = MagicMock()
        mock_file.writelines = self.mock_writelines
        self.mock_open.return_value = open_return
        self.patch2 = patch("os.path.getsize")
        self.mock_size = self.patch2.start()
        self.mock_size.return_value = 0
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        

    def test_can_save_empty_file(self):
        update_data_file("X1", "X")
        self.mock_open.assert_called_with("data/csv/X/X1.csv", "w")
        self.mock_write.assert_called_with("")
    

    def test_can_save_samples(self):
        self.mock_size.return_value = 10
        update_data_file("X1", "X", samples=[{"A": 1, "B": 2}, {"A": 3, "B": 4}])
        self.mock_size.assert_called_with("data/csv/X/X1.csv")
        self.mock_open.assert_called_with("data/csv/X/X1.csv", "a")
        self.mock_writelines.assert_called_with(["1,2\n", "3,4\n"])
        
    
    def test_can_save_samples_with_header(self):
        update_data_file("X1", "X", samples=[{"A": 1, "B": 2}, {"A": 3, "B": 4}])
        self.mock_size.assert_called_with("data/csv/X/X1.csv")
        self.mock_open.assert_called_with("data/csv/X/X1.csv", "a")
        self.mock_writelines.assert_called_with(["A,B\n", "1,2\n", "3,4\n"])



class ResiduesFromModelTests(TestCase):

    def setUp(self):
        self.model = Mock()
        self.residues = [
         {"atomiumId": 1, "atoms": [{"x": 1, "y": 2, "z": 3}]},
         {"atomiumId": 2, "atoms": [{"x": 4, "y": 5, "z": 6}]},
        ]
        self.res1, self.res2, self.res3 = Mock(), Mock(), Mock()
        self.res1.atom.return_value = Mock(location=(1, 2, 3))
        self.res2.atom.return_value = Mock(location=(7, 8, 9))
        self.res3.atom.return_value = Mock(location=(4, 5, 6))
        self.model.residues.side_effect = [
         [self.res1], [self.res2, self.res3]
        ]


    def test_can_get_residues(self):
        residues = get_residues_from_model(self.model, self.residues)
        self.assertEqual(residues, [self.res1, self.res3])



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
        


class ModelToResidueCombinationsTests(TestCase):

    def setUp(self):
        self.patch1 = patch("data.utilities.count_combinations")
        self.mock_count = self.patch1.start()
        self.patch2 = patch("random.sample")
        self.mock_sample = self.patch2.start()
        self.patch3 = patch("data.utilities.split_family")
        self.mock_split = self.patch3.start()
        self.mock_split.side_effect = lambda family: [family]
        self.model = Mock()
        self.model.residues.side_effect = lambda code: []
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()


    def test_can_handle_no_matching_residues(self):
        self.mock_count.return_value = 0
        self.model.residues.side_effect = lambda code: []
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_count.assert_called_with(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, [])
    

    def test_can_handle_insufficient_matching_residues(self):
        self.mock_count.return_value = 0
        self.model.residues.side_effect = lambda code: ["R1", "R2"]
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_count.assert_called_with(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, [])
    

    def test_can_return_single_combination(self):
        self.mock_count.return_value = 1
        residues = Mock(), Mock(), Mock()
        self.model.residues.side_effect = lambda code: residues
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_count.assert_called_with(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, [residues])
    

    def test_can_return_many_combinations(self):
        self.mock_count.return_value = 10
        residues = Mock(), Mock(), Mock(), Mock(), Mock()
        self.model.residues.side_effect = lambda code: residues
        combos = model_to_residue_combos(self.model, "H3")
        self.mock_count.assert_called_with(self.model, "H3")
        self.mock_split.assert_called_with("H3")
        self.model.residues.assert_called_with(code="H")
        self.assertEqual(combos, [
         (residues[0], residues[1], residues[2]),
         (residues[0], residues[1], residues[3]),
         (residues[0], residues[1], residues[4]),
         (residues[0], residues[2], residues[3]),
         (residues[0], residues[2], residues[4]),
         (residues[0], residues[3], residues[4]),
         (residues[1], residues[2], residues[3]),
         (residues[1], residues[2], residues[4]),
         (residues[1], residues[3], residues[4]),
         (residues[2], residues[3], residues[4]),
        ])
    
    def test_can_return_combinations_from_different_subfamilies(self):
        self.mock_count.return_value = 27
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        res = [Mock() for _ in range(9)]
        self.model.residues.side_effect = [res[:3], res[3:6], res[6:]]
        combos = model_to_residue_combos(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, [
         (res[0], res[1], res[3], res[4], res[6], res[7]),
         (res[0], res[1], res[3], res[4], res[6], res[8]),
         (res[0], res[1], res[3], res[4], res[7], res[8]),
         (res[0], res[1], res[3], res[5], res[6], res[7]),
         (res[0], res[1], res[3], res[5], res[6], res[8]),
         (res[0], res[1], res[3], res[5], res[7], res[8]),
         (res[0], res[1], res[4], res[5], res[6], res[7]),
         (res[0], res[1], res[4], res[5], res[6], res[8]),
         (res[0], res[1], res[4], res[5], res[7], res[8]),
         (res[0], res[2], res[3], res[4], res[6], res[7]),
         (res[0], res[2], res[3], res[4], res[6], res[8]),
         (res[0], res[2], res[3], res[4], res[7], res[8]),
         (res[0], res[2], res[3], res[5], res[6], res[7]),
         (res[0], res[2], res[3], res[5], res[6], res[8]),
         (res[0], res[2], res[3], res[5], res[7], res[8]),
         (res[0], res[2], res[4], res[5], res[6], res[7]),
         (res[0], res[2], res[4], res[5], res[6], res[8]),
         (res[0], res[2], res[4], res[5], res[7], res[8]),
         (res[1], res[2], res[3], res[4], res[6], res[7]),
         (res[1], res[2], res[3], res[4], res[6], res[8]),
         (res[1], res[2], res[3], res[4], res[7], res[8]),
         (res[1], res[2], res[3], res[5], res[6], res[7]),
         (res[1], res[2], res[3], res[5], res[6], res[8]),
         (res[1], res[2], res[3], res[5], res[7], res[8]),
         (res[1], res[2], res[4], res[5], res[6], res[7]),
         (res[1], res[2], res[4], res[5], res[6], res[8]),
         (res[1], res[2], res[4], res[5], res[7], res[8]),
        ])
    

    def test_can_handle_insufficient_residues_in_one_subfamily(self):
        self.mock_count.return_value = 0
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        self.model.residues.side_effect = [["H1", "H2", "H3"], ["C1", "C2", "C3"], []]
        combos = model_to_residue_combos(self.model, "H2C2E2")
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, [])
    

    def test_can_limit_number_of_combinations_returned(self):
        self.mock_count.return_value = 27
        self.mock_sample.return_value = [1, 2, 5]
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        res = [Mock() for _ in range(9)]
        self.model.residues.side_effect = [res[:3], res[3:6], res[6:]]
        combos = model_to_residue_combos(self.model, "H2C2E2", limit=3)
        self.mock_count.assert_called_with(self.model, "H2C2E2")
        self.mock_sample.assert_called_with(range(27), 3)
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, [
         (res[0], res[1], res[3], res[4], res[6], res[8]),
         (res[0], res[1], res[3], res[4], res[7], res[8]),
         (res[0], res[1], res[3], res[5], res[7], res[8]),
        ])
    

    def test_can_ignore_seen_ids(self):
        self.mock_count.return_value = 27
        self.mock_sample.return_value = [1, 2, 5]
        self.mock_split.side_effect = lambda family: ["H2", "C2", "E2"]
        res = [Mock() for _ in range(9)]
        self.model.residues.side_effect = [res[:3], res[3:6], res[6:]]
        combos = model_to_residue_combos(self.model, "H2C2E2", limit=3, ignore=[set([res[0].id, res[1].id, res[3].id, res[4].id, res[6].id, res[8].id])])
        self.mock_count.assert_called_with(self.model, "H2C2E2")
        self.mock_sample.assert_called_with(range(27), 3)
        self.mock_split.assert_called_with("H2C2E2")
        self.model.residues.assert_any_call(code="H")
        self.model.residues.assert_any_call(code="C")
        self.model.residues.assert_any_call(code="E")
        self.assertEqual(combos, [
         (res[0], res[1], res[3], res[4], res[7], res[8]),
         (res[0], res[1], res[3], res[5], res[7], res[8]),
        ])



class ModelCombinationsCountTests(TestCase):

    def setUp(self):
        self.patch1 = patch("data.utilities.split_family")
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



class FamilySplittingTests(TestCase):

    def test_can_split_simple_family(self):
        self.assertEqual(split_family("H3"), ["H3"])
        self.assertEqual(split_family("P13"), ["P13"])
    

    def test_can_split_simple_families_with_different_residues(self):
        self.assertEqual(split_family("C4H2"), ["C4", "H2"])
        self.assertEqual(split_family("A12B6C176D1"), ["A12", "B6", "C176", "D1"])



class DatasetSplittingTests(TestCase):

    def test_can_split_dataset(self):
        data = StringIO("ID,col1,col2,positive\nA,1,2,1\nB,10,20,-1\nC,3,4,1")
        df = pd.read_csv(data)
        unlabelled, positives, negatives, core = split_dataset(df)
        self.assertEqual(unlabelled.shape, (3, 3))
        self.assertEqual(positives.shape, (2, 3))
        self.assertEqual(negatives.shape, (1, 3))
        self.assertEqual(core.shape, (3, 2))


