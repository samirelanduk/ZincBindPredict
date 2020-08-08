from atomium import Residue, Atom, Chain, Model, Ligand
from unittest import TestCase
from unittest.mock import patch
from data.common import *

class SequenceSiteToVectorTests(TestCase):

    @patch("data.common.average_hydrophobicity")
    def test_can_get_features_from_sequence(self, mock_hydro):
        mock_hydro.side_effect = [1.1, 0.5, 1.2, 0.8]
        sequence = "xxxxXxxxxxxxxXxx"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8, "hydrophobicity_1": 1.1, "hydrophobicity_3": 0.5
        })
        mock_hydro.assert_any_call(sequence, window=1)
        mock_hydro.assert_any_call(sequence, window=3)
        sequence = "xxxxXxxxxxxxxXxxXXxxxXxxxxxxX"
        self.assertEqual(sequence_site_to_vector(sequence), {
            "gap1": 8, "gap2": 2, "gap3": 0, "gap4": 3, "gap5": 6,
            "hydrophobicity_1": 1.2, "hydrophobicity_3": 0.8
        })
        mock_hydro.assert_any_call(sequence, window=1)
        mock_hydro.assert_any_call(sequence, window=3)



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
    

    def test_average_multiple_hydrophobicity_one_residue(self):
        self.assertEqual(average_hydrophobicity("dxaCwcla", window=3), -0.25)



class StructureSiteToVectorTests(TestCase):

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
        

    def test_can_get_vector_dict(self):
        sample = structure_family_site_to_vector((self.res1, self.res2, self.res3, self.res4))
        self.assertEqual(sample.keys(), {
         "ca_mean", "ca_std", "ca_max", "ca_min", "cb_mean", "cb_std",
         "cb_max", "cb_min", "helix", "strand"
        })
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



class StructureHalfSiteToVectorTests(TestCase):

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
        

    def test_can_get_half_vector_dict(self):
        sample = structure_half_family_site_to_vector((self.res1, self.res2, self.res3, self.res4))
        self.assertEqual(sample.keys(), {
         "ca_mean", "ca_std", "ca_max", "ca_min", "cb_mean", "cb_std",
         "cb_max", "cb_min", "helix", "strand"
        })
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



class ModelToGridTests(TestCase):

    def test_can_convert_model_to_grid(self):
        atom1 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom2 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom3 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom4 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom5 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom6 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        atom7 = Atom("O", 5.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        res1 = Residue(atom1, atom2, atom3)
        res2 = Residue(atom4, atom5, atom6)
        lig = Ligand(atom7)
        model = Model(Chain(res1, res2), lig)
        self.assertEqual(set(model_to_grid(model)), {
            (1, -1, 2), (4, 2, 2), (4, 4, 1), (1, 7, 1),
            (1, 4, 4), (1, 5, 2), (-2, 2, 2), (-2, 4, 1),
            (1, 4, -2), (1, 1, 1), (1, 2, -1), (1, 2, 5)
        })



class LocationToVectorTests(TestCase):

    def test_can_get_location_vector(self):
        atom1 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom2 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom3 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom4 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom5 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom6 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        atom7 = Atom("O", 5.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        atom8 = Atom("ZN", 6.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        res1 = Residue(atom1, atom2, atom3)
        res2 = Residue(atom4, atom5, atom6)
        lig = Ligand(atom7, atom8)
        model = Model(Chain(res1, res2), lig)
        self.assertEqual(location_to_vector((0, 0, 0), model), {
            "8_atom_count": 5, "8_ched_ratio": 0.8, "center_offset": 3.08,
            "16_atom_count": 7, "16_ched_ratio": 0.57
        })
    

    def test_can_get_location_vector_no_atoms(self):
        model = Model()
        self.assertEqual(location_to_vector((0, 0, 0), model), {
            "8_atom_count": 0, "8_ched_ratio": 0, "center_offset": 0,
            "16_atom_count": 0, "16_ched_ratio": 0
        })



class HalfLocationToVectorTests(TestCase):

    def test_can_get_half_location_vector(self):
        atom1 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom2 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom3 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom4 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom5 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom6 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        atom7 = Atom("O", 5.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        atom8 = Atom("ZN", 6.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        atom11 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom12 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom13 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom14 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom15 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom16 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        res1 = Residue(atom1, atom2, atom3)
        res2 = Residue(atom4, atom5, atom6)
        res11 = Residue(atom11, atom12, atom13)
        res12 = Residue(atom14, atom15, atom16)
        lig = Ligand(atom7, atom8)
        model = Model(Chain(res1, res2, id="A"), Chain(res11, res12, id="B"), lig)
        self.assertEqual(half_location_to_vector((0, 0, 0), model, chain="A"), {
            "8_atom_count": 5, "8_ched_ratio": 0.8, "center_offset": 3.08,
            "16_atom_count": 7, "16_ched_ratio": 0.57
        })
    

    def test_can_get_half_location_vector_no_atoms(self):
        model = Model()
        self.assertEqual(half_location_to_vector((0, 0, 0), model, chain="A"), {
            "8_atom_count": 0, "8_ched_ratio": 0, "center_offset": 0,
            "16_atom_count": 0, "16_ched_ratio": 0
        })
        atom1 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom2 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom3 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom4 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom5 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom6 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        res1 = Residue(atom1, atom2, atom3)
        res2 = Residue(atom4, atom5, atom6)
        model = Model(Chain(res1, res2, id="A"))
        self.assertEqual(half_location_to_vector((0, 0, 0), model, chain="B"), {
            "8_atom_count": 0, "8_ched_ratio": 0, "center_offset": 0,
            "16_atom_count": 0, "16_ched_ratio": 0
        })
    

    def test_can_get_half_location_vector_without_chain(self):
        atom1 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom2 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom3 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom4 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom5 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom6 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        atom7 = Atom("O", 5.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        atom8 = Atom("ZN", 6.1, 1.2, 1.3, "O", 0, 0, 0, 0)
        atom11 = Atom("C", 1.1, 1.2, 1.3, "C", 0, 0, 0, 0)
        atom12 = Atom("N", 1.4, 1.5, 1.6, "N", 0, 0, 0, 0)
        atom13 = Atom("O", 1.3, 1.7, 1.8, "O", 0, 0, 0, 0)
        atom14 = Atom("C", 10.1, 2.2, 1.5, "C", 0, 0, 0, 0)
        atom15 = Atom("C", 2.1, 12.2, 1.3, "C", 0, 0, 0, 0)
        atom16 = Atom("S", 1.1, 3.6, 1.3, "S", 0, 0, 0, 0)
        res1 = Residue(atom1, atom2, atom3)
        res2 = Residue(atom4, atom5, atom6)
        res11 = Residue(atom11, atom12, atom13)
        res12 = Residue(atom14, atom15, atom16)
        lig = Ligand(atom7, atom8)
        model = Model(Chain(res1, res2, id="A"), Chain(res11, res12, id="B"), lig)
        self.assertEqual(half_location_to_vector((0, 0, 0), model), {
            "8_atom_count": 9, "8_ched_ratio": 0.78, "center_offset": 2.93,
            "16_atom_count": 13, "16_ched_ratio": 0.54
        })