import os
import unittest
from repeatdesigner import designer as rd

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, 'pdbs')
filename = os.path.join(DATA_DIR,"3atb.pdb")

class TestRepeatBasic(unittest.TestCase):
    def test_pdb(self):
        tpr_des = rd.Design(pdb=filename)
        self.assertEqual(tpr_des.pdb, filename)

    def test_name(self):
        tpr_des = rd.Design(pdb=filename)
        self.assertEqual(tpr_des.name, os.path.join(DATA_DIR,"3atb"))

    def test_sequence(self):
        tpr_des = rd.Design(pdb=filename)
        self.assertEqual(tpr_des.seq, "NSAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAWYNLGNAYYKQGDYDEAIEYYQKALELDPNNAEAKQDLGNAKQKQG")

    def test_targets(self):
        tpr_des = rd.Design(pdb=filename)
        self.assertEqual(tpr_des.targets, range(119))

#    def test_repeats(self):

if __name__ == '__main__':
    unittest.main()
