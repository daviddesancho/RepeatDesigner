"""

This file is part of the RepeatDesigner package


"""

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
        self.targets = targets
        self.compete = compete

    def _read_pdb(self):
        """ 
        Reads PDB structure


        """

