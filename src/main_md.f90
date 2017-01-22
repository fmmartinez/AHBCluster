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

real(8),dimension(:,:),allocatable :: KMatrix, SMatrix

type(Atom),dimension(:),allocatable :: cluster, cluster_initial
type(Forces) :: force, force_initial
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_initial
type(MdData) :: md
type(BasisFunction),dimension(:),allocatable :: phiCov, phiIon, phi
type(BasisFunction),dimension(:),allocatable :: d2pCov, d2pIon, d2p
type(vsl_stream_state) :: stream

nBasisFunCov = 6
nBasisFunIon = 6
nBasisFun = nBasisFunCov + nBasisFunIon

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

allocate(KMatrix(1:nBasisFun,1:nBasisFun))
allocate(SMatrix(1:nBasisFun,1:nBasisFun))

!call read_force_field_file(cluster)
!call initialize_force_field_explicit_H(cluster_initial)
call initialize_force_field_no_H(cluster_initial)

!call read_config_in_XYZ_file(cluster)
!call generate_positions_with_H(cluster_initial)
call generate_positions(cluster_initial)
call write_generated_initial_positions_XYZ(cluster_initial)

!get distances and vectors in all atoms
call get_distances_and_vectors(cluster_initial,atomPairs_initial)
!call get_force_field_pair_parameters_with_H(cluster_initial,atomPairs_initial)
call get_force_field_pair_parameters(cluster_initial,atomPairs_initial)

call initialize_basis_functions_on_each_well(phiCov,phiIon)
phi(1:nBasisFunCov) = phiCov
phi(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = phiIon
call get_overlap_matrix(phi,SMatrix)
call get_double_derivative_basis_functions_on_each_well(d2pCov,d2pIon)
d2p(1:nBasisFunCov) = d2pCov
d2p(nBasisFunCov+1:nBasisFunCov+nBasisFunIon) = d2pIon
call get_kinetic_energy_matrix(phi,d2p,kMatrix)

print '(12f13.6)', KMatrix
print *,''
print '(12f13.9)', SMatrix
stop

!call update_charges_in_complex_and_pairs(cluster_initial,atomPairs_initial)

call get_all_forces(atomPairs_initial,force_initial)

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
