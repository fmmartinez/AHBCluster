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

integer :: i,j,nAtoms,errcode,nBasisFunCov,nBasisFunIon,nBasisFun
integer :: nMapStates,unit1

real(8),dimension(1:3) :: forceCCoM 
real(8),dimension(:),allocatable :: allEigenVal, lambdaVal
real(8),dimension(:,:),allocatable :: HMatrix, allEigenVec
real(8),dimension(:,:),allocatable :: pqp,lql,htemp

type(Atom),dimension(:),allocatable :: cluster, cluster_initial
type(Forces) :: force
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_initial
type(MdData) :: md
type(QuantumStateData) :: quantum, quantum_initial
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

allocate(phiCov(1:nBasisFunCov))
allocate(phiIon(1:nBasisFunIon))
allocate(d2pCov(1:nBasisFunCov))
allocate(d2pIon(1:nBasisFunIon))
allocate(d2p(1:nBasisFun))

allocate(quantum%phi(1:nBasisFun))

allocate(quantum%phiKphi(1:nBasisFun,1:nBasisFun))
allocate(quantum%phiVsphi(1:nBasisFun,1:nBasisFun))
allocate(quantum%hs(1:nMapStates,1:nMapStates))

allocate(quantum%prAHp(1:nBasisFun,1:nBasisFun))

allocate(quantum%SMatrix(1:nBasisFun,1:nBasisFun))
allocate(HMatrix(1:nBasisFun,1:nBasisFun))
allocate(allEigenVec(1:nBasisFun,1:nBasisFun))
allocate(allEigenVal(1:nBasisFun))
allocate(lambdaVal(1:nBasisFun))

allocate(quantum%gridHSolvent(1:nAtoms))

allocate(quantum%pqAp(1:nBasisFun,1:nBasisFun))
allocate(quantum%pqBp(1:nBasisFun,1:nBasisFun))
allocate(quantum%pirp(1:nAtoms))
do i = 1, nAtoms
   allocate(quantum%pirp(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(quantum%pir3p(1:nAtoms))
do i = 1, nAtoms
   allocate(quantum%pir3p(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(quantum%pcr3p(1:nAtoms))
do i = 1, nAtoms
   allocate(quantum%pcr3p(i)%mat(1:nBasisFun,1:nBasisFun))
end do
allocate(quantum%pAHp(1:nBasisFun,1:nBasisFun))
allocate(quantum%pBHp(1:nBasisFun,1:nBasisFun))

allocate(pqp(1:nBasisFun,1:nBasisFun))
allocate(lql(1:nMapStates,1:nMapStates))

allocate(quantum%rm(1:nMapStates))
allocate(quantum%pm(1:nMapStates))
allocate(quantum%mapFactor(1:nMapStates,1:nMapStates))
allocate(quantum%p1(1:nMapStates))
allocate(quantum%q1(1:nMapStates))
allocate(quantum%mapFactor1(1:nMapStates,1:nMapStates))
allocate(quantum%p2(1:nMapStates))
allocate(quantum%q2(1:nMapStates))
allocate(quantum%mapFactor2(1:nMapStates,1:nMapStates))

allocate(quantum%eigenvalues(1:nMapStates))
allocate(quantum%lambda(1:nBasisFun,1:nMapStates))

allocate(quantum%h(1:nMapStates,1:nMapStates))
allocate(htemp(1:nMapStates,1:nMapStates))
allocate(quantum%vas(1:nMapStates,1:nMapStates))
allocate(quantum%vbs(1:nMapStates,1:nMapStates))
allocate(quantum%vhs(1:nMapStates,1:nMapStates))

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
call get_H_grid_Atoms_pos_and_vec(quantum%gridHSolvent,atomPairs_initial)

!call get_force_field_pair_parameters_with_H(cluster_initial,atomPairs_initial)
call get_force_field_pair_parameters(cluster_initial,atomPairs_initial)

call initialize_basis_functions_on_each_well(phiCov,phiIon)
quantum%phi(1:nBasisFunCov) = phiCov
quantum%phi(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = phiIon

call get_overlap_matrix(quantum%phi,quantum%SMatrix)

call get_double_derivative_basis_functions_on_each_well(d2pCov,d2pIon)
d2p(1:nBasisFunCov) = d2pCov
d2p(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = d2pIon

call get_phi_KineticEnergy_phi_matrix(quantum%phi,d2p,quantum%phiKphi)
call get_phi_Vsubsystem_phi_matrix(quantum%phi,atomPairs_initial(1,2)%rij,quantum%phiVsphi)
HMatrix = quantum%phiKphi + quantum%phiVsphi
call get_subsystem_lambdas(HMatrix,quantum%SMatrix,allEigenVec,allEigenVal)

open (newunit=unit1,file='lambdas.log')
do i = 1, nPointsGrid
   do j = 1, nBasisFun
      lambdaVal(j) = sum(allEigenVec(1:nBasisFun,j)*quantum%phi(1:nBasisFun)%gridPointValue(i))
   end do
   write(unit1,'(i3,24f16.4)') i, (lambdaVal(j),j=1,nBasisFun)
end do
close (unit1)

if (nMapStates > 2) then
   quantum%eigenvalues(1:nMapStates) = allEigenVal(1:nMapStates)
   quantum%lambda(1:nBasisFun,1:nMapStates) = allEigenVec(1:nBasisFun,1:nMapStates)
else if (nMapStates == 2) then
   quantum%eigenvalues(1) = allEigenVal(1)
   quantum%eigenvalues(2) = allEigenVal(3)
   quantum%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,1)
   quantum%lambda(1:nBasisFun,2) = allEigenVec(1:nBasisFun,3)
else if (nMapStates == 1) then
   if (md%singleMap == 1) then
      quantum%eigenvalues(1) = allEigenVal(1)
      quantum%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,1)
   else
      quantum%eigenvalues(1) = allEigenVal(3)
      quantum%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,3)
   end if
else
   stop 'error in number of quantum states classically mapped (check nMapStates)'
end if

!call update_charges_in_complex_and_pairs(cluster_initial,atomPairs_initial)
call get_phi_charge_AB_phi_matrix(quantum)
   !previous subroutine called only once, no need to update
call get_phi_inv_r_HS_phi_matrix(quantum)

!h matrix
call get_lambda_h_lambda_matrix(cluster_initial,atomPairs_initial,quantum)
call make_matrix_traceless(quantum%h,quantum%hTraceN,htemp)
quantum%h = htemp

call get_phi_d_VAH_phi_matrix(quantum)
   !previous subroutine called only once, no need to update
call get_phi_d_VBH_phi_matrix(quantum,atomPairs_initial(1,2)%rij)
call get_phi_inv_r3_HS_phi_matrix(quantum)
call get_phi_rc_inv_r3_HS_phi_matrix(atomPairs_initial(1,2)%rij,quantum)

!lastly additional matrix elements for calculating other quantities
call get_phi_rAH_phi_matrix(quantum)

quantum_initial = quantum

do i = 1, md%nTrajectories
   cluster = cluster_initial
   atomPairs = atomPairs_initial
   quantum = quantum_initial
   
   call generate_velocities(cluster,stream,md%initialEqTempInK)
   call remove_CoM_movement(cluster)

   if (md%appMethod == 1) then   
      call do_mapping_variables_sampling(stream,quantum)
      call get_mapFactor(quantum)
      quantum%p1 = 0d0
      quantum%q1 = 0d0
      quantum%mapFactor1 = 0d0
      quantum%p2 = 0d0
      quantum%q2 = 0d0
      quantum%mapFactor2 = 0d0
   else
      call do_coherent_state_variables_sampling(stream,quantum)
      call get_covarFactor(quantum)
      quantum%rm = 0d0
      quantum%pm = 0d0
      quantum%mapFactor = 0d0
   end if
   
   !forces
   if (md%appMethod == 1) then
      call get_all_forces_pbme(cluster_initial,atomPairs_initial,quantum,force,forceCCoM)
   else
      call get_all_forces_fbts(cluster_initial,atomPairs_initial,quantum,force,forceCCoM)
   end if

   print *, 'Trajectory',i,' start'
   
   print *, 'equilibration start'
   if (md%appMethod == 1) then
      if (md%confinement == 1) then
         call run_thermal_equilibration_pbme_confined_cluster(cluster,atomPairs,quantum,force,forceCCoM,md,stream,i)
      else
         call run_thermal_equilibration_pbme(cluster,atomPairs,quantum,force,forceCCoM,md,stream,i)
      end if
   else
      call run_thermal_equilibration_fbts(cluster,atomPairs,quantum,force,forceCCoM,md,stream,i) 
   end if
   print *, 'equilibration end'
   
   !print *, 'production start'
   !call run_nve_pbme(cluster,atomPairs,quantum,force,forceCCoM,md,i)
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
end program clustermd
