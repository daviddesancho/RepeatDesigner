"""

This file is part of the RepeatDesigner package


"""
import numpy as np

import modeller
#from modeller.optimizers import molecular_dynamics, conjugate_gradients
#from modeller.automodel import autosched

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
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')
    return env

def energy(selection):
    """
    Calculates DOPE energy

    Parameters
    ----------
    selection :
        Modeller selection for energy calculation.

    Returns
    -------
    float
        DOPE energy value.

    """
    return selection.assess_dope()

def do_mc(env, mdl, len_mc=10, beta=1.):
    """
    Run MC optimization


    """
    respos = [x.index for x in mdl.residues]
    restyp = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",\
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    while True:
        mdl = modeller.model(env, file=current)
        mdl.write(file='old.pdb')
        rp = random.choice(respos) # randomly select position
        rt = random.choice(restyp) # randomly select residue
        print "\n Iteration : %i, %i, %s"%(n, rp, rt)
        
        # build model for actual mutation
#        try:
#            mdl = mutate_model(mdl, "%s"%rp, rt)
#            mdl.write(file='mutant.pdb')
#            s = selection(mdl)
#            ener_mut = s.assess_dope()
#        
#            # calculate energy interface BA'
#            s = selection()
#            [s.add(mdl.residues["%s:A"%x]) for x in respos]
#            ener_ba = s.assess_dope()
#
#            # build model for competing mutation    
#            mutations = []
#            for r in repeatA:
#                #print r, mdl.residues["%s:A"%r].pdb_name, tpr_residues[r]
#                if mdl.residues["%s:A"%r].pdb_name != tpr_residues[r]:
#                    mutations.append((r, mdl.residues["%s:A"%r].pdb_name))
#                    #print " divergent residue", r, mdl.residues["%s:A"%r].pdb_name, tpr_residues[r]
#            mdl.read(file="initial.pdb")
#            for r, res in mutations:
#                mdl = mutate_model(mdl, "%s"%(r-34), rt)
#            # calculate energy interface AB
#            s = selection()
#            [s.add(mdl.residues["%s:A"%x]) for x in respos0]
#            ener_ab = s.assess_dope()
#            dener = ener_ba - ener_ab
#            ener = w*ener_mut + (1.-w)*dener
#            if ener < ener_prev:
#                print "### ACCEPT ###"
#                ener_prev = ener #[ener_mut, ener_ab, ener_ba]
#                current = 'mutant.pdb'
#                naccept +=1
#                contribs = [ener, ener_mut, ener_ab, ener_ba]
#            else:
#                mdl_prev =  model(env, file='current.pdb')
#                dener = ener - ener_prev
#                if np.exp(-beta*dener) > np.random.random():
#                    print "*** Boltzmann ACCEPT ***"
#                    current = 'mutant.pdb'
#                    naccept +=1
#                    ener_prev = ener 
#                    contribs = [ener, ener_mut, ener_ab, ener_ba]
#                else:
#                    current = 'old.pdb'
#        except OverflowError:
#            current = 'old.pdb'
#
#        print " Current energy %g\n"%ener_prev
#        energy.append(contribs)

        n +=1
        if n >= len_mc:
#            ener_cum.append(energy)
#            mdl = model(env, file=current)
#            mdl.write(file=modelname + "_run%s"%run + ".pdb")
            break

def model_worker(name=None):
    """
    Simple modelling worker function

    Parameters
    ----------
    name : str
        Name of file used as template.

    """
    env = modeller_main()
    mdl = modeller.model(env, file=name)
    s = modeller.selection(mdl)
    ener0 = energy(s) # calculate initial energy
    mdl.write(file='data/initial.pdb') # save initial state
   
    # setup MC 
    beta = 1./50 # inverse temperature
    len_mc = 100 # length of MC run
    n = 0
    naccept = 0
    ener_prev = ener0

    mdl_final = do_mc(env, mdl, initial='data/initial.pdb', len_mc=len_mc, beta=beta)
    return mdl
