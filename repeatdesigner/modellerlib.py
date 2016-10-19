import modeller
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
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')
    return env

def mutate_model(env, name, mdl1, rp, rt):
    """
    Mutates Modeller protein model
    
    Parameters
    ----------
    mdl1 : Modeller model
        The model that will be mutated.
        
    rp : str
        Residue position.
    
    rt : str
    
    Returns
    -------
    """
    chain = "A"
    # Read the original PDB file and copy its sequence to the alignment array:
    ali = modeller.alignment(env)
    ali.append_model(mdl1, atom_files=name, align_codes=name)

    #set up the mutate residue selection segment
    s = modeller.selection(mdl1.chains[chain].residues[int(rp)])

    #perform the mutate residue operation
    s.mutate(residue_type=rt)
    #get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=name)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    #here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    #yes model_copy is the same file as model.  It's a modeller trick.
    mdl2 = modeller.model(env, file=name)

    #required to do a transfer_res_numb
    #ali.append_model(mdl2, atom_files=name, align_codes=name)
    #transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

    #It is usually necessary to write the mutated sequence out and read it in
    #before proceeding, because not all sequence related information about MODEL
    #is changed by this command (e.g., internal coordinates, charges, and atom
    #types and radii are not updated).
    mdl1.write(file=name+rt+"%s"%str(rp)+'.tmp')
    mdl1.read(file=name+rt+"%s"%str(rp)+'.tmp')

    #set up restraints before computing energy
    #we do this a second time because the model has been written out and read in,
    #clearing the previously set restraints
    make_restraints(mdl1, ali)

    #a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms=1
    sched = modeller.automodel.autosched.loop.make_for_model(mdl1)

    #only optimize the selected residue (in first pass, just atoms in selected
    #residue, in second pass, include nonbonded neighboring atoms)
    #set up the mutate residue selection segment
    s = modeller.selection(mdl1.chains[chain].residues[rp])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()
    s.randomize_xyz(deviation=4.0)
    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

    # feels environment (energy computed on pairs that have at least one member
    # in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)
    energy = s.energy()

    #give a proper name
    #mdl1.write(file=name+rt+rp+'.pdb')

    #delete the temporary file
    os.remove(name+rt+"%s"%rp+'.tmp')
    return mdl1


def make_restraints(mdl, aln):
    """Use homologs and dihedral library for dihedral angle restraints """
    rsr = mdl.restraints
    rsr.clear()
    s = modeller.selection(mdl)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

def optimize(atmsel, sched):
    """ Conjugate gradient """
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = modeller.optimizers.conjugate_gradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def refine(atmsel):
    """ Molecular dynamics """
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = modeller.optimizers.molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                            md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False

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
