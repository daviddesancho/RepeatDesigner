"""

This file is part of the RepeatDesigner package


"""
import designlib

class Design(object):
    """
    A class with everything necessary to define a protein repeat design.

    Attributes
    ----------
    pdb         the PDB file name we use as template

    name        the name of the protein

    targets     the regions to model

    compete     the competing regions when we need to introduce negative design

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

class MonteCarlo(object):
    """
    A class for defining and running an MC optimization

    Attributes
    ----------
    len_mc      The length of the optimization

    beta        The inverse temperature for simulated annealing


    """

    def __init__(self, design, nruns=1, len_mc=100, beta=1.):
        """
        Parameters
        ----------
        design : object
            Instance of design class.

        len_mc : int
            Length of optimization

        beta : float
            Inverse temperature.

        """
        self.design = design
        self.nruns = nruns
        self.len_mc = len_mc
        self.beta = beta
    return results
