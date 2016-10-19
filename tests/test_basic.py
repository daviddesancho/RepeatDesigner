import unittest
import coverage
from repeatdesigner import designer as rd

class TestLoad(unittest.TestCase):
    def setUp(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        pass

    def test_name(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.pdb, "pdbs/1vii.pdb")

    def test_targets(self):
        villin_des = rd.Design(pdb="pdbs/1vii.pdb")
        self.assertEqual(villin_des.targets, None)


if __name__ == '__main__':
    unittest.main()
