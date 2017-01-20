module dynamicsroutines
use definitions
implicit none

contains

subroutine run_thermal_equilibration()
implicit none
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
