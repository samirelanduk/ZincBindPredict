from unittest import TestCase
from unittest.mock import patch, MagicMock, Mock
from server.utilities import *

class JobInitializationTests(TestCase):

    @patch("time.time")
    def test_can_initialize_job(self, mock_time):
        mock_time.return_value = 300
        self.assertEqual(initialize_job("ABC"), {
            "id": "300000", "status": "initializing", "protein": "ABC",
            "sites": [], "rejected_sites": []
        })
    

    @patch("time.time")
    def test_can_initialize_job_with_locations(self, mock_time):
        mock_time.return_value = 300
        self.assertEqual(initialize_job("ABC", locations=True), {
            "id": "300000", "status": "initializing", "protein": "ABC",
            "sites": [], "rejected_sites": [], "locations": [], "rejected_locations": []
        })



class JobLocationTests(TestCase):

    def test_can_get_job_location(self):
        self.assertIn(get_job_location("123"), [
            "server/jobs/123.json", "server\\jobs\\123.json"
        ])



class JobSavingTests(TestCase):

    @patch("builtins.open")
    @patch("json.dump")
    def test_can_save_job_without_status(self, mock_dump, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        save_job({"status": "S", "id": "2"})
        mock_dump.assert_called_with({"status": "S", "id": "2"}, mock_file, indent=4)
    

    @patch("builtins.open")
    @patch("json.dump")
    def test_can_save_job_with_status(self, mock_dump, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        save_job({"status": "S", "id": "2"}, status="T")
        mock_dump.assert_called_with({"status": "T", "id": "2"}, mock_file, indent=4)

        

class StructureSavingTests(TestCase):

    @patch("builtins.open")
    def test_can_save_structure(self, mock_open):
        upload = Mock()
        upload.name = "filename.xyz"
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        save_structure_file(upload, "123")
        self.assertTrue(mock_open.call_args_list[0][0][0].endswith("123.xyz"))
        mock_write.assert_called_once_with(upload.read.return_value)
    

    @patch("builtins.open")
    def test_can_save_structure_with_no_file_extension(self, mock_open):
        upload = Mock()
        upload.name = "filename"
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        save_structure_file(upload, "123")
        self.assertTrue(mock_open.call_args_list[0][0][0].endswith("123"))
        mock_write.assert_called_once_with(upload.read.return_value)



class ArgumentParsingTests(TestCase):

    def test_can_parse_args(self):
        sys.argv[1] = "{}"
        self.assertEqual(parse_arguments(), {})



class JobLoadingTests(TestCase):

    @patch("builtins.open")
    @patch("server.utilities.get_job_location")
    @patch("json.load")
    def test_can_load_job(self, mock_load, mock_get, mock_open):
        
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        self.assertEqual(
            load_job("123"), mock_load.return_value
        )
        mock_get.assert_called_with("123")
        mock_open.assert_called_with(mock_get.return_value)
        mock_load.assert_called_with(mock_file)



class FamilySplittingTests(TestCase):

    def test_can_split_single_residue_family(self):
        self.assertEqual(split_family("H3"), [["H", 3]])


    def test_can_split_multi_residue_family(self):
        self.assertEqual(split_family("H3C10"), [["H", 3], ["C", 10]])



class SequenceFamiliesTests(TestCase):

    def test_can_get_available_sequence_families(self):
        self.assertEqual(get_sequence_families(), ["H3", "C4", "C2H2"])



class SequenceToFamilyInputTests(TestCase):

    @patch("server.utilities.split_family")
    def test_can_find_single_residue_family_sites(self, mock_split):
        mock_split.return_value = [["x", 2]]
        inputs = sequence_to_family_inputs("ABXCXDEXFXG", "X2")
        mock_split.assert_called_with("x2")
        self.assertEqual(set(inputs), {
            "abXcXdexfxg", "abXcxdeXfxg", "abXcxdexfXg", "abxcXdeXfxg", "abxcXdexfXg", "abxcxdeXfXg",
        })
    

    @patch("server.utilities.split_family")
    def test_can_find_multi_residue_family_sites(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        inputs = sequence_to_family_inputs("ABXCXDYYEXFXGYYZX", "X2Y3")
        mock_split.assert_called_with("x2y3")
        self.assertEqual(set(inputs), {
            "abXcXdYYexfxgYyzx", "abXcXdYYexfxgyYzx", "abXcXdYyexfxgYYzx", "abXcXdyYexfxgYYzx",
            "abXcxdYYeXfxgYyzx", "abXcxdYYeXfxgyYzx", "abXcxdYyeXfxgYYzx", "abXcxdyYeXfxgYYzx",
            "abXcxdYYexfXgYyzx", "abXcxdYYexfXgyYzx", "abXcxdYyexfXgYYzx", "abXcxdyYexfXgYYzx",
            "abXcxdYYexfxgYyzX", "abXcxdYYexfxgyYzX", "abXcxdYyexfxgYYzX", "abXcxdyYexfxgYYzX",
            "abxcXdYYeXfxgYyzx", "abxcXdYYeXfxgyYzx", "abxcXdYyeXfxgYYzx", "abxcXdyYeXfxgYYzx",
            "abxcXdYYexfXgYyzx", "abxcXdYYexfXgyYzx", "abxcXdYyexfXgYYzx", "abxcXdyYexfXgYYzx",
            "abxcXdYYexfxgYyzX", "abxcXdYYexfxgyYzX", "abxcXdYyexfxgYYzX", "abxcXdyYexfxgYYzX",
            "abxcxdYYeXfXgYyzx", "abxcxdYYeXfXgyYzx", "abxcxdYyeXfXgYYzx", "abxcxdyYeXfXgYYzx",
            "abxcxdYYeXfxgYyzX", "abxcxdYYeXfxgyYzX", "abxcxdYyeXfxgYYzX", "abxcxdyYeXfxgYYzX",
            "abxcxdYYexfXgYyzX", "abxcxdYYexfXgyYzX", "abxcxdYyexfXgYYzX", "abxcxdyYexfXgYYzX",
        })
    

    @patch("server.utilities.split_family")
    def test_no_matches(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        inputs = sequence_to_family_inputs("ABCDEFGH", "X2Y3")
        mock_split.assert_called_with("x2y3")
        self.assertEqual(inputs, [])
    

    @patch("server.utilities.split_family")
    def test_not_enough_matches(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        inputs = sequence_to_family_inputs("ABCDEFGHYYY", "X2Y3")
        self.assertEqual(inputs, [])
        inputs = sequence_to_family_inputs("ABCDEFGHXXXXXYY", "X2Y3")
        self.assertEqual(inputs, [])



class SequenceSiteToVectorTests(TestCase):

    def test_can_convert_sequence_to_vector(self):
        self.assertEqual(sequence_site_to_vector("aaBcc"), [])



class JobModelLoadingTests(TestCase):

    @patch("atomium.open")
    def test_can_load_atomium_model(self, mock_open):
        model = get_model_for_job("filename.xyz")
        self.assertIn(mock_open.call_args_list[0][0][0], [
            "server/jobs/filename.xyz", "server\\jobs\\filename.xyz"
        ])



class StructureFamiliesTests(TestCase):

    def test_can_get_available_structure_families(self):
        self.assertEqual(get_structure_families(), ["H3", "C4", "C2H2"])



class ModelToFamilyInputsTests(TestCase):

    @patch("server.utilities.split_family")
    def test_can_find_single_residue_family_sites(self, mock_split):
        mock_split.return_value = [["x", 2]]
        model = Mock()
        model.residues.return_value = [1, 2, 3]
        inputs = model_to_family_inputs(model, "X2")
        mock_split.assert_called_with("x2")
        model.residues.assert_called_with(code="X")
        self.assertEqual(inputs, [[1, 2], [1, 3], [2, 3]])
    

    @patch("server.utilities.split_family")
    def test_can_find_multi_residue_family_sites(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        model = Mock()
        model.residues.side_effect = [[1, 2, 3], [10, 20, 30, 40]]
        inputs = model_to_family_inputs(model, "X2Y3")
        mock_split.assert_called_with("x2y3")
        model.residues.assert_any_call(code="X")
        model.residues.assert_any_call(code="Y")
        self.assertEqual(inputs, [
            [1, 2, 10, 20, 30], [1, 2, 10, 20, 40],
            [1, 2, 10, 30, 40], [1, 2, 20, 30, 40],
            [1, 3, 10, 20, 30], [1, 3, 10, 20, 40],
            [1, 3, 10, 30, 40], [1, 3, 20, 30, 40],
            [2, 3, 10, 20, 30], [2, 3, 10, 20, 40],
            [2, 3, 10, 30, 40], [2, 3, 20, 30, 40]
        ])
    

    @patch("server.utilities.split_family")
    def test_no_matches(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        model = Mock()
        model.residues.side_effect = [[], []]
        inputs = model_to_family_inputs(model, "X2Y3")
        mock_split.assert_called_with("x2y3")
        self.assertEqual(inputs, [])
    

    @patch("server.utilities.split_family")
    def test_not_enough_matches(self, mock_split):
        mock_split.return_value = [["x", 2], ["y", 3]]
        model = Mock()
        model.residues.side_effect = [[], [10, 20, 30, 40]]
        inputs = model_to_family_inputs(model, "X2Y3")
        self.assertEqual(inputs, [])
        model.residues.side_effect = [[1, 2], [10, 20]]
        inputs = model_to_family_inputs(model, "X2Y3")
        self.assertEqual(inputs, [])



class StructureFamilySiteToVectorTests(TestCase):

    def test_can_convert_structure_family_site_to_vector(self):
        model = Mock()
        self.assertEqual(structure_family_site_to_vector(model), [])



class StructureHalfFamiliesTests(TestCase):

    def test_can_get_available_structure_half_amilies(self):
        self.assertEqual(get_structure_half_families(), ["H2", "C2"])



class StructureFamilyHalfSiteToVectorTests(TestCase):

    def test_can_convert_structure_family_half_site_to_vector(self):
        model = Mock()
        self.assertEqual(structure_family_half_site_to_vector(model), [])



class ModelLocationTests(TestCase):
    
    def test_can_get_model_locations(self):
        model = Mock()
        self.assertEqual(get_model_locations(model), [
            [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
            [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]
        ])



class StructureLocationToVectorTests(TestCase):

    def test_can_convert_location_to_vector(self):
        model = Mock()
        self.assertEqual(structure_location_to_vector([0, 0, 0], model), [])



class StructureLocationToHalfVectorTests(TestCase):

    def test_can_convert_location_to_half_vector(self):
        model = Mock()
        self.assertEqual(structure_location_to_half_vector([0, 0, 0], model), [])