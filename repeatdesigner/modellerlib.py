"""

This file is part of the RepeatDesigner package.

Much of the code has been adapted from Modeller's scripts.
https://salilab.org/modeller/tutorial/basic.html

"""
import os
import numpy as np
import modeller
import modeller.scripts
import modeller.automodel
import modeller.optimizers

def modeller_main():
    """ 
    Basic environment definitions for MODELLER

    """
    modeller.log.none()
    # Set a different value for rand_seed to get a different final model
    env = modeller.environ(rand_seed=-np.random.randint(0,100000))
    env.io.hetatm = True
    #soft sphere potential
    env.edat.dynamic_sphere=False
    env.edat.nonbonded_sel_atoms = 2
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39
    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')
    return env

def gen_align(env, pdb=None, mut=None, out=None):
    """
    Aligns template and test sequences

    Parameters
    ----------
    pdb : str
        PDB filename for template.

    mut : str
        Temporary filename for mutant sequence.

    out : str
        Filename for output.

    Returns
    -------
    aln : object
        Instance of Modeller's alignment class.

    """
    aln = modeller.alignment(env)
    mdl = get_model(env, file=pdb)
    aln.append_model(mdl, align_codes=pdb, atom_files=pdb)
    aln.append(file=mut, align_codes=mut, alignment_format='FASTA')
    aln.align2d()
    aln.write(file=out, alignment_format='PIR')
#    aln.write(file=aliout, alignment_format='PAP')
    return aln

def get_automodel(env, aliin, mut, pdb):
    """
    Generates model from sequence

    env : object
        Modeller environment.

    aliin

    mut

    pdb

    """
    a = modeller.automodel.automodel(env, alnfile=aliin, knowns=pdb, sequence=mut, assess_methods = modeller.automodel.assess.DOPE)
    a.make()
    return a

def get_model(env, file=None):
    """ 
    Generate model using Modeller

    env : object
        Modeller environment.

    """
    return modeller.model(env, file=file)

def get_selection(mdl, sel=None):
    """
    Generate selection

    mdl : object
        Modeller model.

    sel : str
        Selection.

    """
    if not sel:
        return modeller.selection(mdl)
    else:
        s = modeller.selection()
        return modeller.selection(sel)

def get_energy(selection):
    """
    Calculates DOPE score
    TO DO: LJ and coulomb energy

    Parameters
    ----------
    selection : object
        Modeller selection for energy calculation.

    Returns
    -------
    energy : float 
        Modeller energy value.

    """
    #energy = selection.energy(edat=modeller.energy_data(dynamic_sphere=True, \
    #        dynamic_lennard=True, dynamic_coulomb=True, \
    #        nonbonded_sel_atoms=2))
    #return energy[0]
    return selection.assess_dope()

def write_model(mdl, file=None):
    """
    Writes model to file

    Parameters
    ----------
    file : str
        The output file.

    """
    mdl.write(file)

def complete(env, file=None):
    """
    Run complete_pdb[1]_


    Parameters
    ----------
    env

    file : str
        PDB file to read

    Returns
    -------
    mdl :
        Modeller instance of model.

    Notes
    ----------
    ..[1] complete_pdb(): read a PDB file, and fill in any missing residues
    https://salilab.org/modeller/9v1/manual/node343.html
    """
    mdl = modeller.scripts.complete_pdb(env, file)
    return mdl
