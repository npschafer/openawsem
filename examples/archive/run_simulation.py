from openmmtools.integrators import *
from openmmawsem import *

# simulation platform
simulation_platform = "OpenCL" # OpenCL, CUDA, CPU, or Reference

# type of calculation
run_simulation = True
calculate_order_parameters = True
calculate_perturbed_energies = True

# select and prepare input pdb
pdb_id = "1r69"
pdb_chain = 'A'
pdb = "%s.pdb" % pdb_id
# if pdb file is not in the directory, download the pdb file
if not os.path.isfile(pdb):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb.split('.')[0].lower(), pdir='.', file_format='pdb')
    os.rename("pdb%s.ent" % pdb_id, "%s.pdb" % pdb_id)
    
input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, pdb_chain)

# input and output file names
output_pdb_trajectory_file = "output.pdb"
checkpoint_file = "checkpnt.chk"
openmm_formatted_pdb = "%s-openmmawsem.pdb" % pdb_id
cleaned_pdb = "%s-cleaned.pdb" % pdb_id
topology = md.load(openmm_formatted_pdb).topology
num_residues = topology.n_residues

# setup single memory
single_memory = ["%s-cleaned.pdb" % pdb_id, pdb_chain, 1, 1, num_residues, 1]

# order parameter calculation setup
compute_mdtraj = True
compute_total_energy = True
energy_columns = range(4,12) # with column indices starting at 1
order_parameter_file = "order_parameters.txt"
mdtraj_order_parameter_file = "mdtraj_order_parameters.pkl"
rmsd_reference_structure = "%s-openmmawsem.pdb" % pdb_id
qvalue_reference_structure = "%s-cleaned.pdb" % pdb_id
qvalue_reference_chain = pdb_chain

# perturbations
total_energy_column = 15 # specify total energy column (with column indices starting at 1)
perturbation_output_file_name = "order_parameters_with_perturbed_energies.txt"

# simulation setup parameters
reporter_frequency = 1000 # how often to report output
checkpoint_reporter_frequency = 1000000
equilibration_steps = 0 # number of steps to allow equilibration with the qbias potential recording data
num_steps = 10000 # total number of steps in the simulation
step_size = 2*femtoseconds # the step size
q_start = 1.0 # starting q
q_end = 1.0 # ending q
q_increment_size = 0.01 # how often to adjust the q
platform = Platform.getPlatformByName(simulation_platform) 
temperature = 300 # set constant temperature for simulation

# setup Hamiltonian

# choose parameters
# overall scaling
k_awsem=0.5
# con term
k_con=50208
# chain term
k_chain=50208
# chi term
k_chi=251.04
# excl term
k_excl=8368
# rama term
k_rama=8.368
# rama proline term
k_rama_proline=8.368
# amh go term
k_amhgo=1.0
amhgo_min_seq_sep=10
amhgo_contact_threshold=0.8*nanometers
amhgo_chain='A'
amhgo_well_width=0.1
# density dependent AM term
k_am_dd=0.135
am_dd_min_seq_sep=3
am_dd_max_seq_sep=9
density_alpha=1.0
density_normalization=2.75
rho0=1.0
density_min_seq_sep=amhgo_min_seq_sep
am_well_width=0.1
memories=[single_memory]
density_only_from_native_contacts=True
density_pdb_file=cleaned_pdb
density_chain_name=pdb_chain
density_native_contact_min_seq_sep=amhgo_min_seq_sep
density_native_contact_threshold=amhgo_contact_threshold
# qbias
k_qbias=0.0
qbias_min_seq_sep=3
qbias_max_seq_sep=np.inf
qbias_reference_structure="%s-cleaned.pdb" % pdb_id
qbias_reference_chain=pdb_chain

if run_simulation:
    # setup system
    oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=k_awsem) # k_awsem is an overall scaling factor that will affect the relevant temperature scales

    # apply forces
    forces = [
        oa.con_term(k_con=k_con),
        oa.chain_term(k_chain=k_chain),
        oa.chi_term(k_chi=k_chi),
        oa.excl_term(k_excl=k_excl),
        oa.rama_term(k_rama=k_rama),
        oa.rama_proline_term(k_rama_proline=k_rama_proline),
        oa.density_dependent_associative_memory_term(memories, k_am_dd=k_am_dd, am_dd_min_seq_sep=am_dd_min_seq_sep, am_dd_max_seq_sep=am_dd_max_seq_sep, density_alpha=density_alpha, density_normalization=density_normalization, rho0=rho0, density_min_seq_sep=density_min_seq_sep, am_well_width=am_well_width, density_only_from_native_contacts=density_only_from_native_contacts, density_pdb_file=density_pdb_file, density_chain_name=density_chain_name, density_native_contact_min_seq_sep=density_native_contact_min_seq_sep, density_native_contact_threshold=density_native_contact_threshold),
        oa.additive_amhgo_term(cleaned_pdb_filename, amhgo_chain, k_amhgo=k_amhgo, amhgo_min_seq_sep=amhgo_min_seq_sep, amhgo_contact_threshold=amhgo_contact_threshold, amhgo_well_width=amhgo_well_width),
        oa.qbias_term(q_start, qbias_reference_structure, qbias_reference_chain, k_qbias=k_qbias, qbias_min_seq_sep=qbias_min_seq_sep, qbias_max_seq_sep=qbias_max_seq_sep)
    ]
    oa.addForces(forces)

    # start simulation
    collision_rate = 5.0 / picoseconds
    integrator = AndersenVelocityVerletIntegrator(temperature, collision_rate, step_size)
    simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
    simulation.context.setPositions(oa.pdb.positions) # set the initial positions of the atoms
    simulation.context.setVelocitiesToTemperature(temperature) # set the initial velocities of the atoms according to the desired starting temperature
    simulation.minimizeEnergy() # first, minimize the energy to a local minimum to reduce any large forces that might be present
    simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature during simulation

    # if the starting q is not equal to the ending q, change the q periodically throughout the simulation
    if not q_start == q_end:
        steps_per_increment = num_steps/(abs(q_end-q_start)/q_increment_size)
        simulation.reporters.append(PDBReporter(output_pdb_trajectory_file, reporter_frequency)) # output PDBs of simulated structures
        simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency)) # save progress during the simulation
        for q in np.linspace(q_start, q_end, num=int(np.abs((q_start-q_end)/q_increment_size)), endpoint=True):
            simulation.context.setParameter("q0", q)
            simulation.step(int(steps_per_increment))
    # if the starting q is equal to the ending q, do all of the steps at that q and perform equilibration while turning off non-backbone terms
    else:
        simulation.context.setParameter("k_amhgo", 0) # turn off amhgo potential for 1/2 of equilibration
        simulation.context.setParameter("k_am_dd", 0) # turn off am potential for 1/2 of equilibration
        simulation.step(int(0.5*equilibration_steps))
        simulation.context.setParameter("k_amhgo", k_amhgo) # turn on amhgo potential for 1/2 of equilibration
        simulation.context.setParameter("k_am_dd", k_am_dd) # turn on am potential for 1/2 of equilibration   
        simulation.step(int(0.5*equilibration_steps))
        simulation.reporters.append(PDBReporter(output_pdb_trajectory_file, reporter_frequency)) # output PDBs of simulated structures
        simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency)) # save progress during the simulation
        simulation.step(num_steps) # run production simulation

if calculate_order_parameters:
    import pickle
    # when simulation is complete, compute order parameters and energies
    oa = OpenMMAWSEMSystem(openmm_formatted_pdb, k_awsem=k_awsem) # k_awsem is an overall scaling factor that will affect the relevant q scales
    order_parameters_to_compute = {
        "Qvalue": oa.q_value(qvalue_reference_structure, qvalue_reference_chain, min_seq_sep=3, max_seq_sep=np.inf), 
        "LocalQvalue": oa.q_value(qvalue_reference_structure, qvalue_reference_chain, min_seq_sep=3, max_seq_sep=9), 
        "NonlocalQvalue": oa.q_value(qvalue_reference_structure, qvalue_reference_chain, min_seq_sep=10, max_seq_sep=np.inf), 
        "Con": oa.con_term(k_con=k_con),
        "Chain": oa.chain_term(k_chain=k_chain),
        "Chi": oa.chi_term(k_chi=k_chi),
        "Excl": oa.excl_term(k_excl=k_excl),
        "Rama": oa.rama_term(k_rama=k_rama),
        "RamaPro": oa.rama_proline_term(k_rama_proline=k_rama_proline),
        "ddAM": oa.density_dependent_associative_memory_term(memories, k_am_dd=k_am_dd, am_dd_min_seq_sep=am_dd_min_seq_sep, am_dd_max_seq_sep=am_dd_max_seq_sep, density_alpha=density_alpha, density_normalization=density_normalization, rho0=rho0, density_min_seq_sep=density_min_seq_sep, am_well_width=am_well_width, density_only_from_native_contacts=density_only_from_native_contacts, density_pdb_file=density_pdb_file, density_chain_name=density_chain_name, density_native_contact_min_seq_sep=density_native_contact_min_seq_sep, density_native_contact_threshold=density_native_contact_threshold),
        "AMHGo": oa.additive_amhgo_term(cleaned_pdb_filename, amhgo_chain, k_amhgo=k_amhgo, amhgo_min_seq_sep=amhgo_min_seq_sep, amhgo_contact_threshold=amhgo_contact_threshold, amhgo_well_width=amhgo_well_width),
    }

    order_parameters, mdtraj_order_parameters = compute_order_parameters(openmm_formatted_pdb, output_pdb_trajectory_file, order_parameters_to_compute.values(), platform_name=simulation_platform, k_awsem=k_awsem, compute_mdtraj=compute_mdtraj, rmsd_reference_structure=rmsd_reference_structure, compute_total_energy=compute_total_energy, energy_columns=energy_columns)
    order_parameter_names = list(order_parameters_to_compute.keys())
    if compute_mdtraj:
        if rmsd_reference_structure != None:
            order_parameter_names += ["RMSD"]
        order_parameter_names += ["HBondEnergy", "Rg"]
    if compute_total_energy:
        order_parameter_names += ["TotalE"]
    order_parameter_dict = dict(zip(order_parameter_names, range(len(order_parameter_names))))
    np.savetxt(order_parameter_file, order_parameters.T, header=' '.join(order_parameter_dict.keys()), fmt='%.3f')
    pickle.dump(mdtraj_order_parameters, open(mdtraj_order_parameter_file, 'wb'))

if calculate_perturbed_energies:
    oa = OpenMMAWSEMSystem(openmm_formatted_pdb, k_awsem=k_awsem)
    # read order parameters file
    order_parameter_values = np.loadtxt(order_parameter_file)
    # specify term to perturb: column in order parameters file and energy object with all parameters
    perturbations = [
        [(10, oa.density_dependent_associative_memory_term(memories, k_am_dd=k_am_dd, am_dd_min_seq_sep=am_dd_min_seq_sep, am_dd_max_seq_sep=am_dd_max_seq_sep, density_alpha=density_alpha, density_normalization=density_normalization, rho0=rho0, density_min_seq_sep=density_min_seq_sep, am_well_width=am_well_width, density_only_from_native_contacts=density_only_from_native_contacts, density_pdb_file=density_pdb_file, density_chain_name=density_chain_name, density_native_contact_min_seq_sep=density_native_contact_min_seq_sep, density_native_contact_threshold=density_native_contact_threshold)),     
        ],
        [(10, oa.density_dependent_associative_memory_term(memories, k_am_dd=k_am_dd, am_dd_min_seq_sep=am_dd_min_seq_sep, am_dd_max_seq_sep=am_dd_max_seq_sep, density_alpha=density_alpha, density_normalization=density_normalization, rho0=rho0, density_min_seq_sep=density_min_seq_sep, am_well_width=am_well_width, density_only_from_native_contacts=density_only_from_native_contacts, density_pdb_file=density_pdb_file, density_chain_name=density_chain_name, density_native_contact_min_seq_sep=density_native_contact_min_seq_sep, density_native_contact_threshold=density_native_contact_threshold)), 
        ],
    ]
    perturbed_energies = compute_perturbed_energies(openmm_formatted_pdb, output_pdb_trajectory_file, perturbations, order_parameter_values, platform_name=simulation_platform, k_awsem=k_awsem)
    # output new total energies to new order parameters file 
    np.savetxt(perturbation_output_file_name, np.concatenate((order_parameter_values, perturbed_energies), axis=1), fmt='%.3f')
