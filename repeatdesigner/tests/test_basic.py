import os
import unittest
from repeatdesigner import designer as rd
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, 'pdbs')
filename = os.path.join(DATA_DIR,"1vii.pdb")


class TestBasic(unittest.TestCase):
    def test_pdb(self):
        villin_des = rd.Design(pdb=filename)
        self.assertEqual(villin_des.pdb, filename)

    def test_name(self):
        villin_des = rd.Design(pdb=filename)
        self.assertEqual(villin_des.name, os.path.join(DATA_DIR,"1vii"))

    def test_sequence(self):
        villin_des = rd.Design(pdb=filename)
        self.assertEqual(villin_des.seq, "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF")

    def test_targets(self):
        villin_des = rd.Design(filename)
        self.assertEqual(villin_des.targets, range(36))

if __name__ == '__main__':
    unittest.main()
