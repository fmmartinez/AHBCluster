module dynamicsroutines
use definitions
implicit none

contains

subroutine run_nve_dynamics(cluster,atomPairs,force,md)
use ioroutines
use stateevaluation
use energycalculation
implicit none
   
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   
   integer :: i,nAtoms
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   real(8) :: instaTempInK
   
   nAtoms = size(cluster)

   do i = 1, md%prodSteps
      call remove_CoM_movement(cluster)

      call velocity_verlet_int_one_timestep(cluster,atomPairs,force,md)

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_production_trajectory(cluster,i)
      end if
      
      call get_total_potential_energy(atomPairs,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
      totalKinEnergy = get_kinetic_energy(cluster)
      totalEnergy = totalKinEnergy + totalPotEnergy
      
      instaTempInK = get_insta_temperature(totalKinEnergy,nAtoms)
      dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
      totalp = get_total_momentum_magnitude(cluster)
      write(222,'(i10,17f12.6)') i, atomPairs(1,2)%rij, atomPairs(1,3)%rij,&
                                 dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,&
                                 totalPotEnergy,totalKinEnergy,&
                                 totalEnergy, totalp, instaTempInK
   
      if (dcscoms > (nAtoms/(2*0.012d0)**(1/3))) then
         print *, 'production finished at', i,'step due to evaporation'
         exit
      end if
   end do

   print *, 'successful production run'

end subroutine run_nve_dynamics

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
      cluster(j)%pos = cluster(j)%pos + dt*cluster(j)%vel
   end do
   
   !positions changed update positions and vectors, apply constraints and update
   !again
   call get_distances_and_vectors(cluster,atomPairs)
   call do_shake(cluster,atomPairs,mdspecs)
   call get_distances_and_vectors(cluster,atomPairs)

   call update_charges_in_complex_and_pairs(cluster,atomPairs)
   
   !get force
   call get_all_forces(atomPairs,force)

   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do

   call do_rattle(cluster,atomPairs,mdspecs)
end subroutine velocity_verlet_int_one_timestep

subroutine do_shake(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md
   
   integer :: i,j,ai,aj
   real(8) :: esig,esig1,omega2,amti,amtj,gamm,gammi,gammj
   real(8),dimension(1:3) :: tempBondVec
   real(8),dimension(:,:),allocatable :: originalBondVec,at_temp

   allocate(originalBondVec(1:md%nBondConstraints,1:3))
   allocate(at_temp(1:md%nBondConstraints*2,1:3))

   do i = 1, md%nBondConstraints
      ai = 3 + (i-1)*2
      aj = 4 + (i-1)*2
      originalBondVec(i,1:3) = pair(ai,aj)%vectorij

      at_temp(ai,1:3) = at(ai)%pos
      at_temp(aj,1:3) = at(aj)%pos
   end do

   do j = 1, maxShakeCycles
      esig = 0d0
      do i = 1, md%nBondConstraints
         ai = 3 + (i-1)*2
         aj = 4 + (i-1)*2
         tempBondVec = at_temp(aj,1:3) - at_temp(ai,1:3)
         esig1 = abs(sum(tempBondVec**2) - constrainedRS**2)
         esig = max(esig,esig1)
      end do

      if (esig <= toleranceConstraints) then
         exit
      else
         do i = 1, md%nBondConstraints
            ai = 3 + (i-1)*2
            aj = 4 + (i-1)*2
            tempBondVec = pair(ai,aj)%vectorij
            omega2 = constrainedRS**2
            amti = 1d0/at(ai)%mass
            amtj = 1d0/at(aj)%mass

            gamm = (sum(tempBondVec**2)-omega2)/(2d0*(amti+amtj)*sum(originalBondVec(i,1:3)*tempBondVec))
            
            gammi = gamm*amti
            at_temp(ai,1:3) = at_temp(ai,1:3) + originalBondVec(i,1:3)*gammi
            gammj = gamm*amtj
            at_temp(aj,1:3) = at_temp(aj,1:3) - originalBondVec(i,1:3)*gammj
         end do
      end if
   end do

   if (j == maxShakeCycles) then
      print *, 'rattle did not converge'
      stop
   else
      do i = 1, md%nBondConstraints
         ai = 3 + (i-1)*2
         aj = 4 + (i-1)*2

         at(ai)%vel = at(ai)%vel + (at_temp(ai,1:3) - at(ai)%pos)/md%timeStep
         at(aj)%vel = at(aj)%vel + (at_temp(aj,1:3) - at(aj)%pos)/md%timeStep
         
         at(ai)%pos = at_temp(ai,1:3)
         at(aj)%pos = at_temp(aj,1:3)
      end do
   end if
end subroutine do_shake

subroutine do_rattle(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md

   integer :: i,j,ai,aj
   logical :: test
   real(8) :: esig,vvv,wwm,gmma
   
   do j = 1, maxRattleCycles
      esig = 0d0
      do i = 1, md%nBondConstraints
         ai = 3 + (i-1)*2
         aj = 4 + (i-1)*2

         vvv = sum( pair(ai,aj)%vectorij*(at(ai)%vel - at(aj)%vel)  )
         esig = max(esig,abs(vvv))
         wwm = 1d0/at(ai)%mass + 1d0/at(aj)%mass
         gmma = vvv/(wwm*pair(ai,aj)%rij**2)

         at(ai)%vel = at(ai)%vel - gmma/at(ai)%mass*pair(ai,aj)%vectorij
         at(aj)%vel = at(aj)%vel + gmma/at(aj)%mass*pair(ai,aj)%vectorij
      end do

      test = (esig <= toleranceConstraints)

      if (test) then
         exit
      else if (j == maxRattleCycles) then
         print *, 'rattle did not converge'
         stop
      end if
   end do
end subroutine do_rattle

end module dynamicsroutines
