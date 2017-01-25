#!/usr/bin/env python

# # Villin simple design
# Here we attempt to do the simplest possible modelling, which just tries to introduce 
# new side chains at random that improve the total energy. We use the simple villin 
# headpiece subdomain as an example. The more advanced usage designed for repeat 
# proteins are not required. We use a favourite protein as a toy model.

import matplotlib.pyplot as plt
import seaborn as sns
from repeatdesigner import designer as rd

# We will optimize a single residue, which happens to have been mutated experimentally 
# into a histidine. In this case things should work so that we get a rather well 
# converged optimization. 

# We create an instance of the Design class, defining as target residue 25 in the sequence.
villin_des = rd.Design(pdb="pdbs/1vii.pdb", targets=[27])

# Then we create the optimizer, passing arguments like the inverse temperature (`beta`) 
# that will determine the acceptance, the length of the run (`len_mc`) and the number of 
# runs (`nruns`, always think about your number of processors).

mc_villin = rd.Optimizer(villin_des, beta=1e-2, len_mc=100, nruns=50)
mc_villin.run_mc()

fig, ax = plt.subplots()
for k,v in mc_villin.models.iteritems():
    ax.plot(v['score'])
ax.set_ylabel('Energy', fontsize=14)
ax.set_xlabel('MC steps', fontsize=14)
plt.show()

import Bio.PDB
import Bio.Seq
import Bio.SeqUtils
import Bio.pairwise2
import Bio.SeqRecord
import Bio.Align.AlignInfo
for k,v in mc_villin.models.iteritems():
    print "%3i %10.2f %s"%(k, v['score'][-1][0], v['seq'])
#        Bio.Seq.Seq(''.join([Bio.SeqUtils.seq1(x.get_resname()) 
#                    for x in v['model'].get_residues()])))
print

sequences = [Bio.SeqRecord.SeqRecord(x['seq']) for k,x in mc_villin.models.iteritems()] 
align =  Bio.Align.MultipleSeqAlignment(sequences)
summary_align = Bio.Align.AlignInfo.SummaryInfo(align)
print " Consensus sequences:\n -------------------"
print " WT ",villin_des.seq
for t in [0.05, 0.1, 0.2, 0.5, 0.9]:
    print "%.2f"%t, summary_align.dumb_consensus(threshold=(t))
