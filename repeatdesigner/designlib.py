"""

This file is part of the RepeatDesigner package

"""
import os
import sys
import tempfile
import itertools
import random
import numpy as np

import Bio.pairwise2
import Bio.SeqUtils
import Bio.SeqIO
import Bio.SeqRecord

import modellerlib as mdlib

def model_mc_worker(mc_input):
    """
    Simple modelling worker function

    Parameters
    ----------
    mc_input : list
        In which we get the input for the MC optimization.
            design : cls
                Instance of the design class. 
            
            beta : float
                Inverse temperature.
        
            len_mc : int
                length of MC run.
        
    """
    # parse input
    run, design, beta, len_mc = parse_mc_input(mc_input)
    pdb = design.name
    targets = design.targets
    if type(design).__name__ == 'Repeat':
        print " I am a repeat protein!"

    # redirect output by keeping track of sys.stdout
    #nb_stdout = sys.stdout # redirect outputnb_stdout = sys.stdout 
    #sys.stdout = open('data/junk%g.out'%run, 'w') 

    # generate modeller environment
    env = mdlib.modeller_main()

    respos = targets #[x.index for x in mdl.residues]
    restyp = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", \
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", \
            "THR", "TRP", "TYR", "VAL"]

    # generate initial model from PDB file
    mdl = mdlib.get_model(env, file=pdb)
    s = mdlib.get_selection(mdl) 
    seq_prev = design.seq
    ener0 = mdlib.get_energy(s) # calculate initial energy
    contribs = [ener0] #, ener_mut, ener_ab, ener_ba]
    ener_prev = ener0

    n = 0
    ener_mc = [contribs]
    while True:
        rp = random.choice(respos) # randomly select position
        rt = random.choice(restyp) # randomly select residue

        # build model for actual mutation
        print seq_prev
        print rp,
        print rt
        seq, ener = gen_all_models(env, design=design, \
                seq_prev=seq_prev, rp=rp, rt=rt)
                
        # standard acceptance rejection criterion
        seq_prev, ener_prev = boltzmann(ener, \
                ener_prev, seq, seq_prev, beta)

        contribs = [ener_prev]
#        current = 'data/old%g.pdb'%run
        print " Current energy %g\n"%ener_prev
        ener_mc.append(contribs)

        n +=1

        if n >= len_mc:
            # wrap up
            # write mutant sequence to file
            Bio.SeqIO.write(Bio.SeqRecord.SeqRecord(seq_prev, id="data/final%g.fasta"%run), \
                    open("data/final%g.fasta"%run, "w"), "fasta")
            # align sequence and template
            align = mdlib.gen_align(env, pdb=design.pdb, mut="data/final%g.fasta"%run, out="data/align%g.fasta"%run)
            # generate model
            mdl = mdlib.get_automodel(env, "data/align%g.fasta"%run, "data/final%g.fasta"%run, design.pdb)

            ali_tf = tempfile.NamedTemporaryFile(prefix='ali_', suffix='.fasta', \
                delete=False)

            #mdl = mdlib.get_model(env, file=current)
            mdlib.write_model(mdl, file="data/final_run%s"%run + ".pdb")
            #os.remove('data/mutant%g.pdb'%run)
            break
#    sys.stdout = nb_stdout # redirect output    
    return ener_mc

def parse_mc_input(mc_input):
    """
    Parses multiprocessing input

    """
    run = mc_input[0]
    design = mc_input[1][0]
    beta = mc_input[1][1]
    len_mc = mc_input[1][2]
    return run, design, beta, len_mc

def gen_models(env, seq_mut=None, pdb=None):
    """ Generate models based on sequence 

    Parameters
    ----------
    seq_mut : object
        Instance of Seq class.

    pdb : str
        PDB filename for template.

    Returns
    -------
    mdl : object
        Instance of Modeller's model class.

    """
    # write mutant sequence to temporary file
    mut_tf = tempfile.NamedTemporaryFile(prefix='mut_', suffix='.fasta', \
            delete=False)
    Bio.SeqIO.write(Bio.SeqRecord.SeqRecord(seq_mut, id=mut_tf.name), \
            mut_tf.name, "fasta")
    mut_tf.close()

    # align sequence and template
    ali_tf = tempfile.NamedTemporaryFile(prefix='ali_', suffix='.fasta', \
            delete=False)
    align = mdlib.gen_align(env, pdb=pdb, mut=mut_tf.name, out=ali_tf.name)

    # generate model
    mdl = mdlib.get_automodel(env, ali_tf.name, mut_tf.name, pdb)
    
    return mdl

def gen_all_models(env, design=None, seq_prev=None, rp=None, rt=None, weights=[1. ,1.]):
    """ Generate mutated and competing models based on sequence

    Includes modelling the mutant for a given sequence and possible 
    competing selections

    Parameters
    ----------
    env :

    design:

    seq_prev : 

    rp : 

    rt : 

    Returns
    -------
    seq_mut : object
        MutableSeq with the mutated sequence.

    ener : float
        The energy of the mutated model.

    """
    seq_mut, ener_mut = gen_mutated_models(env, design=design, seq_prev=seq_prev, rp=rp, rt=rt)
    if type(design).__name__ == 'Repeat':
        ener_comp = gen_interfaces(env, design=design, seq_mut=seq_mut, seq_prev=seq_prev, rp=rp, rt=rt)

    return seq_mut.toseq(), weights[0]*ener_mut + weights[1]*ener_comp

def gen_mutated_models(env, design=None, seq_prev=None, rp=None, rt=None):
    """ Generate mutated model based on sequence

    Includes modelling the mutant for a given sequence and possible 
    competing selections

    Parameters
    ----------
    env :

    design:

    seq_prev : 

    rp : 

    rt : 

    Returns
    -------
    seq_mut : object
        MutableSeq 

    """
    # generate mutation
    seq_mut = mutator(seq_prev, rp, rt)
    print seq_mut

    # generate model for mutant
    mdl = gen_models(env, seq_mut=seq_mut, pdb=design.pdb)

    # calculate energy for whole mutant
    s = mdlib.get_selection(mdl)
    ener = mdlib.get_energy(s)

    return seq_mut, ener 

def gen_interfaces(env, design=None, seq_mut=None, seq_prev=None, rp=None, rt=None):
    """ 
    Generate competing interfaces
    
    """
    print " Initial repeat:" 
    print [seq_prev[x[0]:x[1]+1] for x in design.repeats]
    # check if mutated residue is in repeat
    try:
        # find position in repeat
        print "Finding mutated residue in repeat"
        mut_rep = [x for x in design.repeats if rp in range(x[0],x[1])][0]
    except IndexError:
        print " IndexError : residue not in repeat"
        return 0.
    ind_mut = design.repeats.index(mut_rep)
    index = rp - mut_rep[0] 
    print " Mutated repeat:", ind_mut, mut_rep
    print [seq_mut[x[0]:x[1]+1] for x in design.repeats]
    print " residue:", index

    # Which interfaces to build
    ireps = range(len(design.repeats))
    interfaces = [y for y in itertools.product(ireps,ireps) if \
            ((y[0]!=y[1]) and (y[1]!=y[0]+1) and ind_mut in y)]
    # build competing interfaces
    ener_interfaces = 0.
    for inter in interfaces:
        # Interfaces are built by threading the repeat sequences into the 1st and 2nd repeat structure
        print " Building competing models", inter
        rep0 = design.repeats[inter[0]]
        rep1 = design.repeats[inter[1]]
        rep_prev = design.repeats[0]
        rep_next = design.repeats[1]
        seq_comp = seq_prev.tomutable()
        # replace original sequence by mutated sequence
        seq_comp[rep_prev[0]:rep_prev[1]+1] = seq_mut[rep0[0]:rep0[1]+1]
        seq_comp[rep_next[0]:rep_next[1]+1] = seq_mut[rep1[0]:rep1[1]+1]
        env.edat.nonbonded_sel_atoms = 2
        mdl = gen_models(env, seq_mut=seq_comp, pdb=design.pdb)
        s0 = mdlib.get_selection(mdl, sel=mdl.residue_range(str(rep_prev[0]), str(rep_prev[1]+1)))
        s1 = mdlib.get_selection(mdl, sel=mdl.residue_range(str(rep_next[0]), str(rep_next[1]+1)))
        s0s1 =mdlib.get_selection(mdl, sel=(mdl.residue_range(str(rep_prev[0]), str(rep_prev[1]+1)), mdl.residue_range(str(rep_next[0]), str(rep_next[1]+1)))) 
        mdlib.get_energy(s0)
        mdlib.get_energy(s1)
        mdlib.get_energy(s0s1)
        #ener_interfaces -= mdlib.get_energy(s0)
        print [seq_comp[x[0]:x[1]+1] for x in design.repeats]
        print
    ener_interfaces = 0.
    return ener_interfaces 

def mutator(seq, rp, rt):
    """
Generates mutation in Seq space

    Parameters
    ----------
    seq : object
        BioPython Seq object to mutate.

    rp : int
        Index indicating position to mutate.

    rt : str
        1 letter code with amino acid type to introduce.

    Returns
    -------
    object
        BioPython Seq object with mutated sequence.

    """
    try:
        mutable_seq = seq.tomutable()
    except AttributeError:
        mutable_seq = seq
    mutable_seq[rp] = Bio.SeqUtils.seq1(rt)
    return mutable_seq

def boltzmann(ener, ener_prev, seq, seq_prev, beta):
    """ Boltzmann acceptance - rejection

    Parameters
    ----------
    ener : float
        New energy.

    ener_prev : float
        Old energy.

    seq : object
        BioPython Seq object with new sequence.

    seq_prev : object
        BioPython Seq object with old sequence.

    beta : float
        Inverse temperature.

    """
    if ener < ener_prev:
        seq_prev = seq
        ener_prev = ener
    else:
        dener = ener - ener_prev
        if np.exp(-beta*dener) > np.random.random():
            seq_prev = seq
            ener_prev = ener 

    return seq_prev, ener_prev
