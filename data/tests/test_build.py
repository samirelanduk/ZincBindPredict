import sys; sys.path.append("data")
from unittest import TestCase
from unittest.mock import patch, Mock, MagicMock
from data.generate_structure_data import main as structures_main


class StructureDatasetBuildingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("builtins.print")
        self.patch2 = patch("data.generate_structure_data.tqdm")
        self.patch3 = patch("data.generate_structure_data.update_data_file")
        self.mock_print = self.patch1.start()
        self.mock_tqdm = self.patch2.start()
        self.mock_data = self.patch3.start()
        self.mock_tqdm.side_effect = lambda l: l


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def check_print_statement(self, fragment):
        for call in self.mock_print.call_args_list:
            if fragment in call[0][0]: break
        else:
            raise ValueError(f"No print call with '{fragment}' in")