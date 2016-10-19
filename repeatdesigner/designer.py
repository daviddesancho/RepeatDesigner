"""

This file is part of the RepeatDesigner package


"""
import designlib

class Design(object):
    """
    A class with everything necessary to define a protein design.

    Attributes
    ----------
    pdb         the PDB file name we use as template

    name        the name of the protein

    targets     the regions to model

    compete     the competing regions when we need to introduce negative design

    """

    def __init__(self, pdb=None, targets=None):
        """
        Parameters
        ----------
        pdb : string
            Name of PDB file to use as template.

        targets : list
            Residues to mutate in order to improve global energy.


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


class Repeat(Design):
    """
    A class with everything necessary to define a protein repeat design

    Attributes
    ----------
    repeats     The regions corresponding to protein repeats

    """
    def __init__(self, pdb=None, targets=None, repeats=None):
        """
        Parameters
        ----------
        compete : list
            Lists of residues for negative design.
        """
        Design.__init__(self, pdb=pdb, targets=targets)
        self.repeats = _parse_repeats(repeats)

    def _parse_repeats(repeats)
        """ 
        Checks whether repeats are same length and sequence.
        
        Parameters
        ----------
        repeats : list
            List of lists containing repeat indexes

        """
        return repeats

class MonteCarlo(object):
    """
    A class for defining and running an MC optimization

    Attributes
    ----------
    len_mc      The length of the optimization

    beta        The inverse temperature for simulated annealing


    """

    def __init__(self, nruns=1, len_mc=10, beta=1.):
        """
        Parameters
        ----------
        len_mc : int
            Length of optimization

        beta : float
            Inverse temperature.

        """
        self.nruns = nruns
        self.len_mc = len_mc
        self.beta = beta


def run_mc(design, mc):
    """
    Parallel MC run generator

    Parameters
    ----------
    design : object
        An instance of the Design class.

    mc : object
        An instance of the MonteCarlo class.

    """

    mdl_final, ener_mc = designlib.model_mc_worker(design.name, beta=mc.beta, len_mc=mc.len_mc)
    return mdl_final, ener_mc
