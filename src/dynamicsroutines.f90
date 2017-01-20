module dynamicsroutines
use definitions
implicit none

contains

subroutine run_thermal_equilibration(cluster,atomPairs,force,md,stream)
use ioroutines
use stateevaluation
use energycalculation
use mkl_vsl_type
use mkl_vsl
implicit none

   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   type(vsl_stream_state),intent(in) :: stream
   
   integer :: i,i_old,try,nAtoms
   real(8) :: tempInK
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   type(Atom),dimension(:),allocatable :: cluster_old
   type(Forces) :: force_old
   type(AtomPairData),dimension(:,:),allocatable :: atomPairs_old
   
   nAtoms = size(cluster)

   allocate(cluster_old(1:nAtoms))
   allocate(atomPairs_old(1:nAtoms,1:nAtoms))

   allocate(force_old%inAtom(1:nAtoms))
   allocate(force_old%atomPair(1:nAtoms,1:nAtoms))
   
   cluster_old = cluster
   atomPairs_old = atomPairs
   force_old = force

   try = 1
   i = 1

   i_old = i

   dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
   
   do while (i <= md%eqSteps)

      tempInK = md%initialEqTempInK + (i/md%eqPhaseSteps)*(md%targetTempInK-md%initialEqTempInK)/md%eqPhases
      
      if (mod(i,md%stepFreqEqSave) == 0) then
         cluster_old = cluster
         atomPairs_old = atomPairs
         force_old = force
         i_old = i
      end if

      if (dcscoms > (nAtoms/(2*0.012d0)**(1/3))) then
         cluster = cluster_old
         atomPairs = atomPairs_old
         force = force_old
         print *, try, 'evaporated and failed at', i
         i = i_old
         try = try + 1
         print *, 'restart', try
         call generate_velocities(cluster,stream,tempInK)
      end if

      call remove_CoM_movement(cluster)
      if (mod(i,(i/md%eqPhaseSteps+1)*md%stepFreqVelRescale) == 0) then
         call do_velocity_rescale(cluster,tempInK)
      end if

      call velocity_verlet_int_one_timestep(cluster,atomPairs,force,md)
      i = i + 1

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_equilibration_trajectory(cluster,i)
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
      if (try == md%maxEqTries) exit
   end do

   if (try == md%maxEqTries) print *, 'stopped equilibration after', try,' tries'
end subroutine run_thermal_equilibration

subroutine velocity_verlet_int_one_timestep(cluster,atomPairs,force,mdspecs)
use stateevaluation
use forcecalculation
implicit none
   
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: mdspecs
   
   integer :: j, nAtoms
   real(8) :: dt,hdt

   nAtoms = size(cluster)

   dt = mdspecs%timeStep
   hdt = mdspecs%halfTimeStep

   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do
   
   do j = 1, nAtoms
      cluster(j)%pos = cluster(j)%pos + dt*cluster(j)%vel
   end do

   !positions changed update positions and vectors
   call get_distances_and_vectors(cluster,atomPairs)

   call update_charges_in_complex_and_pairs(cluster,atomPairs)
   
   !get force
   call get_all_forces(atomPairs,force)

   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do
end subroutine velocity_verlet_int_one_timestep

end module dynamicsroutines
