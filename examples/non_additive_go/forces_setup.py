import os
import sys
from simtk.unit import angstrom
from simtk.unit import kilocalorie_per_mole
from non_aditive_qbias import non_aditive_qbias_term
from openmmawsem import *

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *

def non_aditive_go_term(oa, 
                        q0,
                        reference_pdb_file, 
                        reference_chain_name, 
                        k_qbias=100*kilocalorie_per_mole, 
                        non_aditive_exponent=2.5,
                        qbias_min_seq_sep=3, 
                        qbias_max_seq_sep=np.inf, 
                        contact_threshold=0.8*nanometers, 
                        forceGroup=4):
    k_qbias = k_qbias.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    P=non_aditive_exponent
    
    #Get CA-CB atoms
    c_atoms=[c if c!=-1 else oa.ca[i] for i,c in enumerate(oa.cb)]
    
    #Initialize matrices
    r0_ij=np.zeros([oa.natoms]*2)
    gamma_ij=np.zeros([oa.natoms]*2)
    sigma_ij=np.ones([oa.natoms]*2)*0.1 #No zeros in denominator
    
    #Read native distances
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, 
        reference_pdb_file, reference_chain_name=reference_chain_name,
        min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, 
        contact_threshold=contact_threshold, Qflag=0)
        
    #Set parameters:
    for i,j,(g,r0,s) in structure_interactions:
        if r0>contact_threshold:
            continue
        gamma_ij[i,j]=g
        r0_ij[i,j]=r0.value_in_unit(nanometer)
        sigma_ij[i,j]=s
    
    alpha=gamma_ij.sum(axis=0)
    alpha[alpha<1]=1.0 #No zeros in denominator
    print(f"TOTAL CONTACTS={gamma_ij.sum()}")
    print(f"TEST={(gamma_ij.sum(axis=0)**2).sum()}")

    #Initialize Force
    go = CustomGBForce()
    
    #Add Parameters
    go.addTabulatedFunction("gamma_ij", Discrete2DFunction(oa.natoms,oa.natoms,gamma_ij.T.flatten()))
    go.addTabulatedFunction("r0_ij", Discrete2DFunction(oa.natoms,oa.natoms,r0_ij.T.flatten()))
    go.addTabulatedFunction("sigma_ij", Discrete2DFunction(oa.natoms,oa.natoms,sigma_ij.T.flatten()))
    go.addPerParticleParameter("id",)
    go.addPerParticleParameter("alpha",)
    
    #Compute non-additive go term
    go.addComputedValue("Eij","gamma_ij(id1,id2)*exp(-(r-r0_ij(id1,id2))^2/(2*sigma_ij(id1,id2)^2))",CustomGBForce.ParticlePair)
    go.addComputedValue("Ei",f"1/alpha^{P}*Eij",CustomGBForce.SingleParticle)
    go.addEnergyTerm(f"-0.5*{k_qbias}*Ei^{P}",CustomGBForce.SingleParticle)
    
    #Finish force setting
    go.setForceGroup(forceGroup)
    go.setCutoffDistance(2)
    go.setNonbondedMethod(CustomGBForce.CutoffNonPeriodic)
   
    #Add particles
    for i,a in enumerate(alpha):
        go.addParticle([i,a])
    
    return go
   
def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        non_aditive_go_term(oa, 
                       q0=1,
                       reference_pdb_file=f"1r69-cleaned.pdb",
                       reference_chain_name="A"),
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        contact_term(oa),
        # for membrane protein simulation use contact_term below.
        # contact_term(oa, z_dependent=True, inMembrane=True, membrane_center=membrane_center, k_relative_mem=3),
        beta_term_1(oa),
        beta_term_2(oa),
        beta_term_3(oa),
        pap_term_1(oa),
        pap_term_2(oa),
        # er_term(oa),
        # membrane_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center),
        # membrane_preassigned_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center, zimFile="PredictedZim"),
        # fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
        # debye_huckel_term(oa, chargeFile="charge.txt"),
    ]
    if computeQ:
        forces.append(rg_term(oa))
        forces.append(q_value(oa, "crystal_structure-cleaned.pdb", forceGroup=1))
        # forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", residueIndexGroup=list(range(0, 15)), forceGroup=1))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces
