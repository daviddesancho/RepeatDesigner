"""

This file is part of the RepeatDesigner package


"""
import numpy as np

from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched

def modeller_main():
    """ 
    Basic environment definitions for MODELLER

    """
    log.none()
    # Set a different value for rand_seed to get a different final model
    env = environ(rand_seed=-np.random.randint(0,100000))
    env.io.hetatm = True
    #soft sphere potential
    env.edat.dynamic_sphere=False
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')
    return env

def model_worker(name=None):
    """
    Simple modelling worker function

    Parameters
    ----------
    name : str
        Name of file used as template.

    """
    env = modeller_main()
    mdl = model(env, file=name)
