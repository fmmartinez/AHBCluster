program clustermd
use definitions
use mkl_vsl_type
use mkl_vsl
use ioroutines
use dynamicsroutines
use stateevaluation
use forcecalculation
use quantumcalculations
use maproutines
implicit none

integer :: i,nAtoms,errcode,nBasisFunCov,nBasisFunIon,nBasisFun
integer :: nMapStates

real(8),dimension(1:3) :: forceCCoM, forceCCoM_initial
real(8),dimension(:,:),allocatable :: KMatrix, VhMatrix, SMatrix, HMatrix
real(8),dimension(:,:),allocatable :: pqp,lql,htemp

type(Atom),dimension(:),allocatable :: cluster, cluster_initial
type(Forces) :: force, force_initial
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_initial
type(MdData) :: md
type(QuantumStateData) :: pbme
type(BasisFunction),dimension(:),allocatable :: phiCov, phiIon
type(BasisFunction),dimension(:),allocatable :: d2pCov, d2pIon, d2p
type(vsl_stream_state) :: stream


call read_md_input_file(nAtoms,nBasisFunCov,nBasisFunIon,nMapStates,md)
errcode = vslnewstream(stream,brng,md%seed)

nBasisFun = nBasisFunCov + nBasisFunIon

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
allocate(d2pCov(1:nBasisFunCov))
allocate(d2pIon(1:nBasisFunIon))
allocate(d2p(1:nBasisFun))

allocate(pbme%phi(1:nBasisFun))

allocate(KMatrix(1:nBasisFun,1:nBasisFun))
allocate(VhMatrix(1:nBasisFun,1:nBasisFun))
allocate(SMatrix(1:nBasisFun,1:nBasisFun))
allocate(HMatrix(1:nBasisFun,1:nBasisFun))

allocate(pbme%gridHSolvent(1:nAtoms))

allocate(pbme%pqAp(1:nBasisFun,1:nBasisFun))
allocate(pbme%pqBp(1:nBasisFun,1:nBasisFun))
allocate(pbme%pirp(1:nAtoms))
do i = 1, nAtoms
   allocate(pbme%pirp(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(pbme%pir2p(1:nAtoms,1:3))
do i = 1, nAtoms
   allocate(pbme%pir2p(i,1)%mat(1:nBasisFun,1:nBasisFun))
   allocate(pbme%pir2p(i,2)%mat(1:nBasisFun,1:nBasisFun))
   allocate(pbme%pir2p(i,3)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(pbme%pir3p(1:nAtoms))
do i = 1, nAtoms
   allocate(pbme%pir3p(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(pbme%pcr3p(1:nAtoms))
do i = 1, nAtoms
   allocate(pbme%pcr3p(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(pbme%pAHp(1:nBasisFun,1:nBasisFun))
allocate(pbme%pBHp(1:nBasisFun,1:nBasisFun))

allocate(pqp(1:nBasisFun,1:nBasisFun))
allocate(lql(1:nMapStates,1:nMapStates))

allocate(pbme%rm(1:nMapStates))
allocate(pbme%pm(1:nMapStates))
allocate(pbme%mapFactor(1:nMapStates,1:nMapStates))

allocate(pbme%eigenvalues(1:nMapStates))
allocate(pbme%lambda(1:nBasisFun,1:nMapStates))

allocate(pbme%h(1:nMapStates,1:nMapStates))
allocate(htemp(1:nMapStates,1:nMapStates))
allocate(pbme%vas(1:nMapStates,1:nMapStates))
allocate(pbme%vbs(1:nMapStates,1:nMapStates))
allocate(pbme%vhs(1:nMapStates,1:nMapStates))

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
call get_H_grid_Atoms_pos_and_vec(pbme%gridHSolvent,atomPairs_initial)

!call get_force_field_pair_parameters_with_H(cluster_initial,atomPairs_initial)
call get_force_field_pair_parameters(cluster_initial,atomPairs_initial)

call initialize_basis_functions_on_each_well(phiCov,phiIon)
pbme%phi(1:nBasisFunCov) = phiCov
pbme%phi(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = phiIon
call get_overlap_matrix(pbme%phi,SMatrix)
call get_double_derivative_basis_functions_on_each_well(d2pCov,d2pIon)
d2p(1:nBasisFunCov) = d2pCov
d2p(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = d2pIon
call get_phi_KineticEnergy_phi_matrix(pbme%phi,d2p,KMatrix)
call get_phi_Vsubsystem_phi_matrix(pbme%phi,atomPairs_initial(1,2)%rij,VhMatrix)

HMatrix = KMatrix + VhMatrix
call get_subsystem_lambdas(HMatrix,SMatrix,pbme%lambda,pbme%eigenvalues)

!call update_charges_in_complex_and_pairs(cluster_initial,atomPairs_initial)
call get_phi_charge_AB_phi_matrix(pbme)
   !previous subroutine called only once, no need to update
call get_phi_inv_r_HS_phi_matrix(pbme)
!h matrix
call get_lambda_h_lambda_matrix(cluster_initial,atomPairs_initial,pbme)
call make_matrix_traceless(pbme%h,pbme%hTraceN,htemp)
pbme%h = htemp

call get_phi_d_VAH_phi_matrix(pbme)
   !previous subroutine called only once, no need to update
call get_phi_d_VBH_phi_matrix(pbme,atomPairs_initial(1,2)%rij)
call get_phi_inv_r2_HS_phi_matrix(pbme)
call get_phi_inv_r3_HS_phi_matrix(pbme)
call get_phi_rc_inv_r3_HS_phi_matrix(atomPairs_initial(1,2)%rij,pbme)

!print *, 'a'
!print '(12f9.4)', pbme%pAHp
!print *, 'b'
!print '(12f9.4)', pbme%pBHp
!print *, 'l1'
!print '(12f9.4)', pbme%lambda(1:12,1)
!print *, 'l2'
!print '(12f9.4)', pbme%lambda(1:12,2)

pbme%rm(1) = 0.1205d0
pbme%pm(1) = 0.1244d0
do i = 2, nMapStates
pbme%rm(i) = 0.0d0
pbme%pm(i) = 0.0d0
end do
call get_mapFactor(pbme)
!forces
call get_all_forces_pbme(cluster_initial,atomPairs_initial,pbme,force_initial,forceCCoM_initial)

!call get_phi_q_phi_matrix(phi,gridHSolvent,pqp)
!call get_lambda_q_lambda_matrix(lambda,pqp,lql)
!print *, 'qm',get_map_contribution(lql,mapFactor)

do i = 1, md%nTrajectories
   cluster = cluster_initial
   atomPairs = atomPairs_initial
   force = force_initial
   forceCCoM = forceCCoM_initial
   
   call generate_velocities(cluster,stream,md%initialEqTempInK)
   call remove_CoM_movement(cluster)

   print *, 'Trajectory',i,' start'
   
   print *, 'equilibration start'
   !call run_thermal_equilibration(cluster,atomPairs,force,md,stream,i)
   call run_thermal_equilibration_pbme_only_11(cluster,atomPairs,pbme,force,forceCCoM,md,stream,i)
   print *, 'equilibration end'
   
   !print *, 'production start'
   !call run_nve_dynamics(cluster,atomPairs,force,md,i)
   !print *, 'production end'
   
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
