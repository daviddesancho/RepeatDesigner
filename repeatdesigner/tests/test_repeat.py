import unittest
from repeatdesigner import designer as rd

class TestRepeatBasic(unittest.TestCase):
    def test_pdb(self):
        tpr_des = rd.Design(pdb="repeatdesigner/tests/pdbs/3atb.pdb")
        self.assertEqual(tpr_des.pdb, "repeatdesigner/tests/pdbs/3atb.pdb")

    def test_name(self):
        tpr_des = rd.Design(pdb="repeatdesigner/tests/pdbs/3atb.pdb")
        self.assertEqual(tpr_des.name, "repeatdesigner/tests/pdbs/3atb")

    def test_sequence(self):
        tpr_des = rd.Design(pdb="repeatdesigner/tests/pdbs/3atb.pdb")
        self.assertEqual(tpr_des.seq, "NSAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAKQDLGNAKQKQG")

    def test_targets(self):
        tpr_des = rd.Design(pdb="repeatdesigner/tests/pdbs/3atb.pdb")
        self.assertEqual(tpr_des.targets, range(119))

if __name__ == '__main__':
    unittest.main()
