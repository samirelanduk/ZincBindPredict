from unittest import TestCase
from unittest.mock import patch
from predict.utilities import *

class CategoryFromArgumentsAndFilesystemTests(TestCase):

    @patch("os.listdir")
    def test_can_get_categories_from_file_system(self, mock_list):
        mock_list.side_effect = [
            ["category1", "category2"], ["dataset1.csv", "dataset2.csv"], ["dataset3.csv", "xxx"]
        ]
        categories = get_categories_from_arguments_and_filesystem([])
        self.assertEqual(categories, {
            "category1": {"dataset1": {}, "dataset2": {}},
            "category2": {"dataset3": {}}
        })
        mock_list.assert_any_call("data/csv")
        mock_list.assert_any_call("data/csv/category1")
        mock_list.assert_any_call("data/csv/category2")
    

    @patch("os.listdir")
    def test_can_get_categories_from_file_system_filtered(self, mock_list):
        mock_list.side_effect = [
            ["category1", "category2"], ["dataset1.csv", "dataset2.csv"], ["dataset3.csv", "xxx"]
        ]
        categories = get_categories_from_arguments_and_filesystem([
            "xxx", "--categories=category1", "--datasets=dataset2"
        ])
        self.assertEqual(categories, {
            "category1": {"dataset2": {}}
        })
        mock_list.assert_any_call("data/csv")
        mock_list.assert_any_call("data/csv/category1")
        mock_list.assert_any_call("data/csv/category2")