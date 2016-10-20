import unittest
from repeatdesigner import designer as rd

class TestBasic(unittest.TestCase):
    def test_pdb(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.pdb, "pdbs/1vii.pdb")

    def test_name(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.name, "pdbs/1vii")

    def test_sequence(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.seq, "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF")

    def test_targets(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.targets, range(36))

class TestRepeatBasic(unittest.TestCase):
    def test_pdb(self):
        tpr_des = rd.Design(pdb="pdbs/3atb.pdb")
        self.assertEqual(tpr_des.pdb, "pdbs/3atb.pdb")

    def test_name(self):
        tpr_des = rd.Design(pdb="pdbs/3atb.pdb")
        self.assertEqual(tpr_des.name, "pdbs/3atb")

    def test_sequence(self):
        tpr_des = rd.Design(pdb="pdbs/3atb.pdb")
        self.assertEqual(tpr_des.seq, "NSAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAKQDLGNAKQKQG")
    def test_targets(self):
        tpr_des = rd.Design(pdb="pdbs/3atb.pdb")
        self.assertEqual(tpr_des.targets, range(119))

if __name__ == '__main__':
    unittest.main()
