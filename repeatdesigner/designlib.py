"""

This file is part of the RepeatDesigner package


"""
import os
import sys
import numpy as np
import random
import modellerlib as mdlib

def parse_mc_input(mc_input):
    """
    Parses multiprocessing input

    """
    run = mc_input[0]
    pdb = mc_input[1][0]
    beta = mc_input[1][1]
    len_mc = mc_input[1][2]
    targets = mc_input[1][3]
    return run, pdb, beta, len_mc, targets

def model_mc_worker(mc_input):
    """
    Simple modelling worker function

    Parameters
    ----------
    pdb : str
        Name of file used as template.
    
    beta : float
        Inverse temperature.

    len_mc : int
        length of MC run.

    targets : list
        Residues to mutate

    """
    
    # parse input
    run, pdb, beta, len_mc, targets = parse_mc_input(mc_input)
    nb_stdout = sys.stdout # redirect outputnb_stdout = sys.stdout 
    sys.stdout = open('data/junk%g.out'%run, 'w') # redirect output
    env = mdlib.modeller_main()

    mdl = mdlib.get_model(env, file=pdb)

    respos = targets #[x.index for x in mdl.residues]
    restyp = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",\
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    n = 0
    naccept = 0
    s = mdlib.get_selection(mdl) 
    ener0 = mdlib.get_energy(s) # calculate initial energy
    ener_prev = ener0
    ener_mc = []
    current = pdb
    while True:
        mdl = mdlib.get_model(env, file=current)
        mdlib.write_model(mdl, file='data/old%g.pdb'%run)
        rp = random.choice(respos) # randomly select position
        rt = random.choice(restyp) # randomly select residue
        print "\n Iteration : %i, %i, %s"%(n, rp, rt)
        
        # build model for actual mutation
        try:
            mdl = mdlib.mutate_model(env, pdb, "data/mut%g_"%run , mdl, rp-1, rt)
            mdlib.write_model(mdl, file='data/mutant%g.pdb'%run)
            s = mdlib.get_selection(mdl)
            ener_mut = mdlib.get_energy(s)

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
            ener = ener_mut
            if ener < ener_prev:
                print "### ACCEPT ###"
#                ener_prev = ener #[ener_mut, ener_ab, ener_ba]
                current = 'data/mutant%g.pdb'%run
                naccept +=1
                contribs = [ener] #, ener_mut, ener_ab, ener_ba]
                ener_prev = ener
            else:
                dener = ener - ener_prev
                if np.exp(-beta*dener) > np.random.random():
                    print "*** Boltzmann ACCEPT ***"
                    current = 'data/mutant%g.pdb'%run
                    naccept +=1
                    ener_prev = ener 
                    contribs = [ener] #, ener_mut, ener_ab, ener_ba]
                else:
                    print "### ACCEPT ###"
                    current = 'data/old%g.pdb'%run
                    contribs = [ener_prev]
        except OverflowError:
            current = 'data/old%g.pdb'%run
#
        print " Current energy %g\n"%ener_prev
        ener_mc.append(contribs)

        n +=1
        if n >= len_mc:
#            ener_cum.append(energy)
            mdl = mdlib.get_model(env, file=current)
            mdlib.write_model(mdl, file="data/final_run%s"%run + ".pdb")
            os.remove('data/old%g.pdb'%run)
            os.remove('data/mutant%g.pdb'%run)
            break
    sys.stdout = nb_stdout # redirect output    
    return ener_mc
