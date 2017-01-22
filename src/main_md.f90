program clustermd
use definitions
use mkl_vsl_type
use mkl_vsl
use ioroutines
use dynamicsroutines
use stateevaluation
use energycalculation
use forcecalculation
implicit none


integer :: i,i_old,j,unit1,nAtoms,seed,errcode,try
integer :: eqSteps,stepFreqEqSave,stepFreqOutTrajectory
integer :: maxEqTries

real(8) :: f1,f2,f3,f4,fr,rah,rab,rbh,totalEnergy,dcscoms
real(8) :: totalp,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy,totalKinEnergy
real(8) :: tempInK,tempInK_old
real(8),dimension(1:3) :: uvah,uvab,uvhb,fa,fb,fh

type(Atom),dimension(:),allocatable :: cluster,cluster_old
type(Forces) :: force, force_old
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_old
type(MdData) :: md
type(vsl_stream_state) :: stream

call read_md_input_file(nAtoms,md)
errcode = vslnewstream(stream,brng,md%seed)

allocate(cluster(1:nAtoms),cluster_old(1:nAtoms))
allocate(atomPairs(1:nAtoms,1:nAtoms),atomPairs_old(1:nAtoms,1:nAtoms))

allocate(force%inAtom(1:nAtoms),force_old%inAtom(1:nAtoms))
allocate(force%atomPair(1:nAtoms,1:nAtoms),force_old%atomPair(1:nAtoms,1:nAtoms))

!call read_force_field_file(cluster)
call initialize_force_field_explicit_H(cluster)
!call read_config_in_XYZ_file(cluster)
call generate_positions(cluster)
call write_generated_initial_positions_XYZ(cluster)

call generate_velocities(cluster,stream,md%initialEqTempInK)
call remove_CoM_movement(cluster)

!get distances and vectors in all atoms
call get_distances_and_vectors(cluster,atomPairs)
call get_force_field_pair_parameters(cluster,atomPairs)

call update_charges_in_complex_and_pairs(cluster,atomPairs)

call get_all_forces(atomPairs,force)

print *, 'equilibration start'
call run_thermal_equilibration(cluster,atomPairs,force,md,stream)
print *, 'equilibration end'
print *, 'production start'
call run_nve_dynamics(cluster,atomPairs,force,md)
print *, 'production end'

end program clustermd
