from unittest import TestCase
from unittest.mock import patch
from server.utilities import sequence_to_family_inputs

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