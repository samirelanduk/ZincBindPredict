from unittest.mock import patch, call
from django.test import TestCase
from server.sequence_job import main as seq_main

class SequenceJobTests(TestCase):

    def setUp(self):
        self.patch1 = patch("server.sequence_job.parse_arguments")
        self.mock_parse = self.patch1.start()
        self.mock_parse.return_value = {"sequence": "abc", "job_id": "12"}
        self.patch2 = patch("server.sequence_job.load_job")
        self.mock_load = self.patch2.start()
        self.mock_load.return_value = {"id": "1"}
        self.patch3 = patch("server.sequence_job.get_sequence_families")
        self.mock_fam = self.patch3.start()
        self.mock_fam.return_value = ["X1", "Y1"]
        self.patch4 = patch("server.sequence_job.save_job")
        self.mock_save = self.patch4.start()
        self.patch5 = patch("server.sequence_job.sequence_to_family_inputs")
        self.mock_inputs = self.patch5.start()
        self.patch5.return_value = ["aaa", "bbb"]
        self.patch6 = patch("server.sequence_job.sequence_site_to_vector")
        self.mock_vectors = self.patch6.start()
        self.patch6.side_effect = lambda i: i.upper()
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()
        self.patch4.stop()
        self.patch5.stop()
        self.patch6.stop()


    def test_successful_job(self):
        seq_main()
        self.mock_load.assert_called_with("12")
        self.assertEqual(
            self.mock_save.call_args_list[0],
            call({"id": "1"}, status="Looking for X1 sites")
        )
        self.assertEqual(self.mock_inputs.call_args_list[0], call("ABC", "X1"))