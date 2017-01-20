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

real(8) :: f1,f2,f3,f4,fr,rah,rab,rbh,totalEnergy,dcscoms
real(8) :: totalp,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy,totalKinEnergy
real(8) :: tempInK,tempInK_old
real(8),dimension(1:3) :: uvah,uvab,uvhb,fa,fb,fh

type(Atom),dimension(:),allocatable :: cluster,cluster_old
type(Forces) :: force, force_old
type(AtomPairData),dimension(:,:),allocatable :: atomPairs, atomPairs_old
type(MdData) :: mdspecs
type(vsl_stream_state) :: stream

seed = 2
errcode = vslnewstream(stream,brng,seed)

nAtoms = 203 !cluster 3, solvent 9x2
mdspecs%timeStep = 0.001d0
mdspecs%halfTimeStep = mdspecs%timeStep/2d0
tempInK = 10

allocate(cluster(1:nAtoms),cluster_old(1:nAtoms))
allocate(atomPairs(1:nAtoms,1:nAtoms),atomPairs_old(1:nAtoms,1:nAtoms))

allocate(force%inAtom(1:nAtoms),force_old%inAtom(1:nAtoms))
allocate(force%atomPair(1:nAtoms,1:nAtoms),force_old%atomPair(1:nAtoms,1:nAtoms))

!call read_force_field_file(cluster)
call initialize_force_field_explicit_H(cluster)
!call read_config_in_XYZ_file(cluster)
call generate_positions(cluster)
call generate_velocities(cluster,stream,tempInK)
call remove_CoM_movement(cluster)

!get distances and vectors in all atoms
call get_distances_and_vectors(cluster,atomPairs)
call get_force_field_pair_parameters(cluster,atomPairs)

call update_charges_in_complex_and_pairs(cluster,atomPairs)

call get_all_forces(atomPairs,force)

dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)

open(newunit=unit1,file='ini_generated.xyz')
write(unit1,*) nAtoms
write(unit1,*) 'symbol - positions x y z --',i
do j = 1, nAtoms
   write(unit1,'(a4,3f14.8)') cluster(j)%symbol, cluster(j)%pos(1:3)
end do
close(unit1)

cluster_old = cluster
atomPairs_old = atomPairs
force_old = force
tempInK_old = tempInK

try = 1
i = 1

i_old = i

do while (i <= 375000)

   tempInK = 10d0 + (i/37500)*14
   
   if (mod(i,37500) == 0) then
      cluster_old = cluster
      atomPairs_old = atomPairs
      force_old = force
      i_old = i
   end if

   if (dcscoms > 7d0) then
      cluster = cluster_old
      atomPairs = atomPairs_old
      force = force_old
      print *, try, 'failed at', i
      i = i_old
      try = try + 1
      print *, 'restart', try
      call generate_velocities(cluster,stream,tempInK)
   end if

   if (mod(i,(i/37500+1)*10) == 0) then
      call remove_CoM_movement(cluster)
      call do_velocity_rescale(cluster,tempInK)
   end if

   call velocity_verlet_int_one_timestep(cluster,atomPairs,force,mdspecs)
   i = i + 1

   if (mod(i,250) == 0) then
      open(newunit=unit1,file='trajectory.xyz',position='append')
      write(unit1,*) nAtoms
      write(unit1,*) 'symbol - positions x y z --',i
      do j = 1, nAtoms
         write(unit1,'(a4,3f14.8)') cluster(j)%symbol, cluster(j)%pos(1:3)
      end do
      close(unit1)
   end if

   call get_total_potential_energy(atomPairs,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
   totalKinEnergy = get_kinetic_energy(cluster)
   totalEnergy = totalKinEnergy + totalPotEnergy
   
   dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
   totalp = get_total_momentum_magnitude(cluster)
   write(111,'(i10,16f12.6)') i, atomPairs(1,2)%rij, atomPairs(1,3)%rij,&
                              dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,&
                              totalPotEnergy,totalKinEnergy,&
                              totalEnergy, totalp
   if (try == 10) exit
end do

if (try == 10) print *, 'stopped after 10 tries'

end program clustermd
