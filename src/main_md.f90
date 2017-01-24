program clustermd
use definitions
use mkl_vsl_type
use mkl_vsl
use ioroutines
use dynamicsroutines
use stateevaluation
use forcecalculation
use quantumcalculations
implicit none

integer :: i,nAtoms,errcode,nBasisFunCov,nBasisFunIon,nBasisFun
integer :: nMapStates

real(8),dimension(:),allocatable :: eigenvalues
real(8),dimension(:,:),allocatable :: KMatrix, VhMatrix, SMatrix, HMatrix
real(8),dimension(:,:),allocatable :: pqAp,pqBp
real(8),dimension(:,:),allocatable :: lambda,h

type(Atom),dimension(:),allocatable :: cluster, cluster_initial
type(Forces) :: force, force_initial
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_initial
type(MdData) :: md
type(BasisFunction),dimension(:),allocatable :: phiCov, phiIon, phi
type(BasisFunction),dimension(:),allocatable :: d2pCov, d2pIon, d2p
type(EvalOnGridHData),dimension(:),allocatable :: gridHSolvent
type(MatrixList),dimension(:),allocatable :: pirp
type(vsl_stream_state) :: stream

nBasisFunCov = 6
nBasisFunIon = 6
nBasisFun = nBasisFunCov + nBasisFunIon

nMapStates = 3

call read_md_input_file(nAtoms,md)
errcode = vslnewstream(stream,brng,md%seed)

allocate(cluster(1:nAtoms))
allocate(atomPairs(1:nAtoms,1:nAtoms))
allocate(cluster_initial(1:nAtoms))
allocate(atomPairs_initial(1:nAtoms,1:nAtoms))

allocate(force%inAtom(1:nAtoms))
allocate(force%atomPair(1:nAtoms,1:nAtoms))
allocate(force_initial%inAtom(1:nAtoms))
allocate(force_initial%atomPair(1:nAtoms,1:nAtoms))

allocate(phiCov(1:nBasisFunCov))
allocate(phiIon(1:nBasisFunIon))
allocate(phi(1:nBasisFun))
allocate(d2pCov(1:nBasisFunCov))
allocate(d2pIon(1:nBasisFunIon))
allocate(d2p(1:nBasisFun))

allocate(gridHSolvent(1:nAtoms))

allocate(KMatrix(1:nBasisFun,1:nBasisFun))
allocate(VhMatrix(1:nBasisFun,1:nBasisFun))
allocate(SMatrix(1:nBasisFun,1:nBasisFun))
allocate(HMatrix(1:nBasisFun,1:nBasisFun))

allocate(pqAp(1:nBasisFun,1:nBasisFun))
allocate(pqBp(1:nBasisFun,1:nBasisFun))
allocate(pirp(1:nAtoms))
do i = 1, nAtoms
   allocate(pirp(i)%mat(1:nBasisFun,1:nBasisFun))
end do

allocate(eigenvalues(1:nMapStates))
allocate(lambda(1:nBasisFun,1:nMapStates))
allocate(h(1:nMapStates,1:nMapStates))

!call read_force_field_file(cluster)
!call initialize_force_field_explicit_H(cluster_initial)
call initialize_force_field_no_H(cluster_initial)

!call read_config_in_XYZ_file(cluster)
!call generate_positions_with_H(cluster_initial)
call generate_positions(cluster_initial)
call write_generated_initial_positions_XYZ(cluster_initial)

!get distances and vectors in all atoms
call get_distances_and_vectors(cluster_initial,atomPairs_initial)
!get distances and vectors of all solvent atoms with H grid
call get_H_grid_Atoms_pos_and_vec(gridHSolvent,atomPairs_initial)

!call get_force_field_pair_parameters_with_H(cluster_initial,atomPairs_initial)
call get_force_field_pair_parameters(cluster_initial,atomPairs_initial)

call initialize_basis_functions_on_each_well(phiCov,phiIon)
phi(1:nBasisFunCov) = phiCov
phi(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = phiIon
call get_overlap_matrix(phi,SMatrix)
call get_double_derivative_basis_functions_on_each_well(d2pCov,d2pIon)
d2p(1:nBasisFunCov) = d2pCov
d2p(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = d2pIon
call get_phi_KineticEnergy_phi_matrix(phi,d2p,KMatrix)
call get_phi_Vsubsystem_phi_matrix(phi,atomPairs_initial(1,2)%rij,VhMatrix)

HMatrix = KMatrix + VhMatrix
call get_subsystem_lambdas(HMatrix,SMatrix,lambda,eigenvalues)

!call update_charges_in_complex_and_pairs(cluster_initial,atomPairs_initial)
call get_phi_q_ABS_phi_matrix(phi,pqAp,pqBp)
call get_phi_inv_r_HS_phi_matrix(phi,gridHSolvent,pirp)

call get_lambda_h_lambda_matrix(cluster_initial,atomPairs_initial,&
                                 eigenvalues,lambda,pqAp,pqBp,pirp,h)
print '(3f13.4)', eigenvalues
print *,''
print *, h
stop

!call get_all_forces(atomPairs_initial,force_initial)

do i = 1, md%nTrajectories
   cluster = cluster_initial
   atomPairs = atomPairs_initial
   force = force_initial
   
   call generate_velocities(cluster,stream,md%initialEqTempInK)
   call remove_CoM_movement(cluster)

   print *, 'Trajectory',i,' start'
   
   print *, 'equilibration start'
   call run_thermal_equilibration(cluster,atomPairs,force,md,stream,i)
   print *, 'equilibration end'
   
   print *, 'production start'
   call run_nve_dynamics(cluster,atomPairs,force,md,i)
   print *, 'production end'
   
   print *, 'Trajectory',i,' end'
   print *, ' '
end do

deallocate(cluster)
deallocate(atomPairs)
deallocate(cluster_initial)
deallocate(atomPairs_initial)

deallocate(force%inAtom)
deallocate(force%atomPair)
deallocate(force_initial%inAtom)
deallocate(force_initial%atomPair)
end program clustermd
