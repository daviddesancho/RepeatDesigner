"""

This file is part of the RepeatDesigner package


"""
import designlib

class Design(object):
    """
    A class with everything necessary to build a designed repeat
    protein.

    Attributes
    ----------


    """
    def __init__(self, pdb=None, targets=None, compete=None):
        """
        Parameters
        ----------
        pdb : string
            Name of PDB file to use as template.

        targets : list
            Residues to mutate in order to improve global energy.

        compete : list
            Lists of residues for negative design.

        """
        self.pdb = pdb
        self.name = self._pdb_name()
        self.targets = targets
        self.compete = compete

    def _pdb_name(self):
        ie = self.pdb.rfind('.pdb')
        return self.pdb[:ie]


    def _parse_targets(self):
        """ 
        Parse target regions for modelling


        """
        if self.targets is None:
            pass

    def launch_model(self):
        """
        Launches external optimization of defined models.

        """
        results = designlib.model_worker(name=self.name)
        return results
