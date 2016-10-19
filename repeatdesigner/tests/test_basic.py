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
        self.assertEqual(villin_des.targets, None)

    def test_targets_high(self):



if __name__ == '__main__':
    unittest.main()
