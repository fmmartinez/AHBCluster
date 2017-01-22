program clustermd
use definitions
use mkl_vsl_type
use mkl_vsl
use ioroutines
use dynamicsroutines
use stateevaluation
use forcecalculation
implicit none

integer :: i,nAtoms,errcode

type(Atom),dimension(:),allocatable :: cluster, cluster_initial
type(Forces) :: force, force_initial
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_initial
type(MdData) :: md
type(vsl_stream_state) :: stream

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
