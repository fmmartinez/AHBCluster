module dynamicsroutines
use definitions
implicit none

contains

subroutine run_nve_dynamics(cluster,atomPairs,force,md,trj)
use ioroutines
use stateevaluation
use energycalculation
implicit none
   integer,intent(in) :: trj
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   
   character(17) :: outLogFile, trjxyzFile
   integer :: i,nAtoms,unit1,unit2
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   real(8) :: instaTempInK
   
   nAtoms = size(cluster)
   
   write(outLogFile,'(a9,i4.4,a4)') 'pd-output',trj,'.log'
   write(trjxyzFile,'(a9,i4.4,a4)') 'prod-traj',trj,'.xyz'
   open(newunit=unit1,file=outLogFile)
   open(newunit=unit2,file=trjxyzFile)
   
   do i = 1, md%prodSteps
      if (mod(i,md%stepFreqComRemoval) == 0 ) call remove_CoM_movement(cluster)

      call velocity_verlet_int_one_timestep(cluster,atomPairs,force,md)

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_xyz_trajectory(cluster,i,unit2)
      end if
      
      if (mod(i,md%stepFreqOutLog) == 0) then      
         call get_total_potential_energy(atomPairs,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
         totalKinEnergy = get_kinetic_energy(cluster)
         totalEnergy = totalKinEnergy + totalPotEnergy
         
         instaTempInK = get_insta_temperature(totalKinEnergy,nAtoms,md%nBondConstraints)
         dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
         totalp = get_total_momentum_magnitude(cluster)
         write(unit1,'(i10,17f12.6)') i, atomPairs(1,2)%rij, atomPairs(1,3)%rij,&
                                    dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,&
                                    totalPotEnergy,totalKinEnergy,&
                                    totalEnergy, totalp, instaTempInK
      end if

      if (maxval(atomPairs(1,1:nAtoms)%rij) > (2d0*(nAtoms/(0.012d0))**(1d0/3d0))) then
         print *, 'production finished at', i,'step due to evaporation',&
                  maxval(atomPairs(1,1:nAtoms)%rij)
         exit
      end if
   end do

   if (i >= md%prodSteps) then
      print *, 'finished successful production run', i
   else
      print *, 'production run finished early', i
   end if
   
   close(unit2)
   close(unit1)
end subroutine run_nve_dynamics

subroutine run_thermal_equilibration(cluster,atomPairs,force,md,stream,trj)
use ioroutines
use stateevaluation
use energycalculation
use mkl_vsl_type
use mkl_vsl
implicit none
   integer,intent(in) :: trj
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   type(vsl_stream_state),intent(in) :: stream
   
   character(17) :: outLogFile, trjxyzFile
   integer :: i,i_old,try,nAtoms,unit1,unit2
   real(8) :: tempInK, maxDistAS, clusterRadius
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
   clusterRadius = (nAtoms/(0.012d0))**(1d0/3d0)
   
   write(outLogFile,'(a9,i4.4,a4)') 'eq-output',trj,'.log'
   write(trjxyzFile,'(a9,i4.4,a4)') 'equi-traj',trj,'.xyz'
   open(newunit=unit1,file=outLogFile)
   open(newunit=unit2,file=trjxyzFile)
   
   do while (i <= md%eqSteps)

      tempInK = md%initialEqTempInK + (i/md%eqPhaseSteps)*(md%targetTempInK-md%initialEqTempInK)/md%eqPhases
      
      if (mod(i,md%stepFreqEqSave) == 0) then
         cluster_old = cluster
         atomPairs_old = atomPairs
         force_old = force
         i_old = i
      end if
      
      maxDistAS = clusterRadius*1.75d0*(0.75d0+(i/md%eqPhaseSteps)*(0.25/md%eqPhases))
      if (maxval(atomPairs(1,1:nAtoms)%rij) > maxDistAS) then
         cluster = cluster_old
         atomPairs = atomPairs_old
         force = force_old
         print *, try, 'evaporated and failed at', i, maxval(atomPairs(1,1:nAtoms)%rij)
         i = i_old
         try = try + 1
         print *, 'restart', try
         call generate_velocities(cluster,stream,tempInK)
         call remove_CoM_movement(cluster)
         call do_rattle(cluster,atomPairs,md)
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
      end if

      if (mod(i,md%stepFreqCoMremoval) == 0 ) call remove_CoM_movement(cluster)
      if (mod(i,(i/md%eqPhaseSteps+1)*md%stepFreqVelRescale) == 0) then
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
         call remove_CoM_movement(cluster)
         call do_rattle(cluster,atomPairs,md)
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
      end if

      call velocity_verlet_int_one_timestep(cluster,atomPairs,force,md)
      i = i + 1

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_xyz_trajectory(cluster,i,unit2)
      end if

      if (mod(i,md%stepFreqOutLog) == 0) then      
         call get_total_potential_energy(atomPairs,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
         totalKinEnergy = get_kinetic_energy(cluster)
         totalEnergy = totalKinEnergy + totalPotEnergy
         
         dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
         totalp = get_total_momentum_magnitude(cluster)
         write(unit1,'(i10,16f12.6)') i, atomPairs(1,2)%rij, atomPairs(1,3)%rij,&
                                    dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,&
                                    totalPotEnergy,totalKinEnergy,&
                                    totalEnergy, totalp
      end if

      if (try > md%maxEqTries) exit
   end do

   if (try > md%maxEqTries) then
      print *, 'stopped equilibration after', try,' tries'
      stop
   end if

   close(unit2)
   close(unit1)
end subroutine run_thermal_equilibration

subroutine run_nve_pbme(cluster,atomPairs,p,force,forceCCoM,md,trj)
use ioroutines
use stateevaluation
use energycalculation
use maproutines
implicit none
   integer,intent(in) :: trj
   real(8),dimension(1:3),intent(inout) :: forceCCoM
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   type(QuantumStateData),intent(inout) :: p
   
   character(17) :: outLogFile1, outLogFile2, trjxyzFile1
   integer :: i,j,try,nAtoms,nMapStates,unit1,unit2,unit3
   real(8) :: tempInK, maxDistAS, clusterRadius
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   real(8) :: solPol, instaTempInK
   
   nAtoms = size(cluster)
   nMapStates = size(p%rm)

   clusterRadius = (nAtoms/(0.012d0))**(1d0/3d0)
   
   write(outLogFile1,'(a9,i4.4,a4)') 'outputPd1',trj,'.log'
   write(outLogFile2,'(a9,i4.4,a4)') 'outputPd2',trj,'.log'
   write(trjxyzFile1,'(a9,i4.4,a4)') 'xtrajProd',trj,'.xyz'
   open(newunit=unit1,file=outLogFile1)
   open(newunit=unit2,file=outLogFile2)
   open(newunit=unit3,file=trjxyzFile1)
   
   do i = 1, md%prodSteps
      if (mod(i,md%stepFreqCoMremoval) == 0 ) call remove_CoM_movement(cluster)
      
      call velocity_verlet_int_one_timestep_pbme(cluster,atomPairs,p,force,forceCCoM,md)

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_xyz_trajectory(cluster,i,unit3)
      end if

      if (mod(i,md%stepFreqOutLog) == 0) then
         if (atomPairs(1,2)%rij /= atomPairs(1,2)%rij) stop 'NaN found'      
         call get_total_potential_energy_pbme(atomPairs,p,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
         totalKinEnergy = get_kinetic_energy(cluster)
         totalEnergy = totalKinEnergy + totalPotEnergy
         
         instaTempInK = get_insta_temperature(totalKinEnergy,nAtoms,md%nBondConstraints)
         dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
         totalp = get_total_momentum_magnitude(cluster)
         solPol = get_solvent_polarization(cluster,atomPairs)
         write(unit1,'(i10,24f12.6)') i, atomPairs(1,2)%rij, solpol,&
            dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess, totalPotEnergy,totalKinEnergy,&
            totalEnergy, totalp, instaTempInK
         write(unit2,'(i10,24f12.6)') i, (p%h(j,j),j=1,nMapStates),&
            ((p%rm(j)**2+p%pm(j)**2)/(2d0*hbar),j=1,nMapStates)
      end if
      
      if (maxval(atomPairs(1,1:nAtoms)%rij) > 2d0*clusterRadius) then
         print *, 'production finished at', i,'step due to evaporation',&
                  maxval(atomPairs(1,1:nAtoms)%rij)
         exit
      end if
   end do

   if (i >= md%prodSteps) then
      print *, 'finished successful production run', i
   else
      print *, 'production run finished early', i
   end if

   close(unit3)
   close(unit2)
   close(unit1)
end subroutine run_nve_pbme

subroutine run_thermal_equilibration_pbme(cluster,atomPairs,p,force,forceCCoM,md,stream,trj)
use ioroutines
use stateevaluation
use energycalculation
use maproutines
use mkl_vsl_type
use mkl_vsl
implicit none
   integer,intent(in) :: trj
   real(8),dimension(1:3),intent(inout) :: forceCCoM
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   type(QuantumStateData),intent(inout) :: p
   type(vsl_stream_state),intent(in) :: stream
   
   character(17) :: outLogFile1, outLogFile2, trjxyzFile1
   integer :: i,i_old,j,k,try,nAtoms,nMapStates,unit1,unit2,unit3
   real(8) :: tempInK, distCS, clusterRadius
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   real(8) :: solPol
   real(8),dimension(1:3) :: forceCCoM_old, com
   type(Atom),dimension(:),allocatable :: cluster_old
   type(Forces) :: force_old
   type(QuantumStateData) :: p_old
   type(AtomPairData),dimension(:,:),allocatable :: atomPairs_old
   
   nAtoms = size(cluster)
   nMapStates = size(p%rm)

   allocate(cluster_old(1:nAtoms))
   allocate(atomPairs_old(1:nAtoms,1:nAtoms))

   allocate(force_old%inAtom(1:nAtoms))
   allocate(force_old%atomPair(1:nAtoms,1:nAtoms))
   
   cluster_old = cluster
   atomPairs_old = atomPairs
   force_old = force
   p_old = p

   try = 1
   i = 1

   i_old = i

   dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
   clusterRadius = 2d0*(nAtoms/0.012d0)**(1d0/3d0)
   
   write(outLogFile1,'(a9,i4.4,a4)') 'outputEq1',trj,'.log'
   write(outLogFile2,'(a9,i4.4,a4)') 'outputEq2',trj,'.log'
   write(trjxyzFile1,'(a9,i4.4,a4)') 'xtrajEqui',trj,'.xyz'
   open(newunit=unit1,file=outLogFile1)
   open(newunit=unit2,file=outLogFile2)
   open(newunit=unit3,file=trjxyzFile1)
   
   do while (i <= md%eqSteps)

      tempInK = md%initialEqTempInK + (i/md%eqPhaseSteps)*(md%targetTempInK-md%initialEqTempInK)/md%eqPhases
      
      if (mod(i,md%stepFreqEqSave) == 0) then
         cluster_old = cluster
         atomPairs_old = atomPairs
         force_old = force
         forceCCoM_old = forceCCoM
         i_old = i
         p_old = p
      end if
      
      call get_center_of_mass_vector(cluster,com)
      do k = 1, nAtoms
         distCS = sqrt(sum((com - cluster(k)%pos)**2))
         if (distCS > clusterRadius) then
            cluster = cluster_old
            atomPairs = atomPairs_old
            force = force_old
            forceCCoM_old = forceCCoM
            print *, try, 'evaporated and failed at', i, distCS
            i = i_old
            p = p_old
            try = try + 1
            print *, 'restart', try
            call generate_velocities(cluster,stream,tempInK)
            call remove_CoM_movement(cluster)
            call do_rattle(cluster,atomPairs,md)
            call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
            !regenerate regenerable stuff from quantum state
            !p%rm = 0d0
            !p%pm = 0d0
            !p%rm(1) = 0.1205d0
            !p%pm(1) = 0.1244d0
            call do_mapping_variables_sampling(stream,p)
            call get_mapFactor(p)
         end if
      end do

      if (mod(i,md%stepFreqCoMremoval) == 0 ) call remove_CoM_movement(cluster)
      if (mod(i,(i/md%eqPhaseSteps+1)*md%stepFreqVelRescale) == 0) then
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
         call remove_CoM_movement(cluster)
         call do_rattle(cluster,atomPairs,md)
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
      end if

      call velocity_verlet_int_one_timestep_pbme(cluster,atomPairs,p,force,forceCCoM,md)
      i = i + 1

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_xyz_trajectory(cluster,i,unit3)
      end if

      if (mod(i,md%stepFreqOutLog) == 0) then
         if (atomPairs(1,2)%rij /= atomPairs(1,2)%rij) stop 'NaN found'      
         call get_total_potential_energy_pbme(atomPairs,p,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
         totalKinEnergy = get_kinetic_energy(cluster)
         totalEnergy = totalKinEnergy + totalPotEnergy
         
         dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
         totalp = get_total_momentum_magnitude(cluster)
         solPol = get_solvent_polarization(cluster,atomPairs)
         write(unit1,'(i10,15f12.6)') i, atomPairs(1,2)%rij, solpol,&
            dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess, totalPotEnergy,totalKinEnergy,&
            totalEnergy, totalp
         write(unit2,'(i10,24f12.6)') i, (p%h(j,j),j=1,nMapStates),&
            ((p%rm(j)**2+p%pm(j)**2)/(2d0*hbar),j=1,nMapStates)
         !write(999,'(i10,24f12.6)') i, (force%inAtom(j)%total,j=1,nAtoms)
         !write(888,'(i10,24f12.6)') i, (atomPairs(3,j)%rij,j=1,nAtoms)
         write(887,*) i, get_apparent_rAH(p)
      end if

      if (try > md%maxEqTries) exit
   end do

   if (try > md%maxEqTries) then
      print *, 'stopped equilibration after', try,' tries'
      stop
   end if

   close(unit3)
   close(unit2)
   close(unit1)
end subroutine run_thermal_equilibration_pbme

subroutine run_thermal_equilibration_pbme_confined_cluster(cluster,atomPairs,p,force,forceCCoM,md,stream,trj)
use ioroutines
use stateevaluation
use energycalculation
use maproutines
use mkl_vsl_type
use mkl_vsl
implicit none
   integer,intent(in) :: trj
   real(8),dimension(1:3),intent(inout) :: forceCCoM
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: md
   type(QuantumStateData),intent(inout) :: p
   type(vsl_stream_state),intent(in) :: stream
   
   character(17) :: outLogFile1, outLogFile2, trjxyzFile1
   integer :: i,j,try,nAtoms,nMapStates,unit1,unit2,unit3
   real(8) :: tempInK, maxDistAS, clusterRadius, clusterMaxRadius
   real(8) :: dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess
   real(8) :: totalPotEnergy,totalKinEnergy,totalEnergy, totalp
   real(8) :: solPol
   
   nAtoms = size(cluster)
   nMapStates = size(p%rm)

   dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
   
   write(outLogFile1,'(a9,i4.4,a4)') 'outputEq1',trj,'.log'
   write(outLogFile2,'(a9,i4.4,a4)') 'outputEq2',trj,'.log'
   write(trjxyzFile1,'(a9,i4.4,a4)') 'xtrajEqui',trj,'.xyz'
   open(newunit=unit1,file=outLogFile1)
   open(newunit=unit2,file=outLogFile2)
   open(newunit=unit3,file=trjxyzFile1)
   
   do i = 1, md%eqSteps
      tempInK = md%initialEqTempInK + (i/md%eqPhaseSteps)*(md%targetTempInK-md%initialEqTempInK)/md%eqPhases

      if (mod(i,md%stepFreqCoMremoval) == 0 ) call remove_CoM_movement(cluster)
      if (mod(i,(i/md%eqPhaseSteps+1)*md%stepFreqVelRescale) == 0) then
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
         call remove_CoM_movement(cluster)
         call do_rattle(cluster,atomPairs,md)
         call do_velocity_rescale(cluster,tempInK,md%nBondConstraints)
      end if

      call velocity_verlet_int_one_timestep_pbme_confined_cluster(cluster,atomPairs,p,force,forceCCoM,md)

      if (mod(i,md%stepFreqOutTrajectory) == 0) then
         call write_xyz_trajectory(cluster,i,unit3)
      end if

      if (mod(i,md%stepFreqOutLog) == 0) then
         if (atomPairs(1,2)%rij /= atomPairs(1,2)%rij) stop 'NaN found'      
         call get_total_potential_energy_pbme(atomPairs,p,ec,ecslj,ecsel,ecs,esslj,essel,essb,ess,totalPotEnergy)
         totalKinEnergy = get_kinetic_energy(cluster)
         totalEnergy = totalKinEnergy + totalPotEnergy
         
         dcscoms = get_distance_solvent_CoM_complex_CoM(cluster)
         totalp = get_total_momentum_magnitude(cluster)
         solPol = get_solvent_polarization(cluster,atomPairs)
         write(unit1,'(i10,15f12.6)') i, atomPairs(1,2)%rij, solpol,&
            dcscoms, ec,ecslj,ecsel,ecs,esslj,essel,essb,ess, totalPotEnergy,totalKinEnergy,&
            totalEnergy, totalp
         write(unit2,'(i10,24f12.6)') i, (p%h(j,j),j=1,nMapStates),&
            ((p%rm(j)**2+p%pm(j)**2)/(2d0*hbar),j=1,nMapStates)
      end if
   end do

   close(unit3)
   close(unit2)
   close(unit1)
end subroutine run_thermal_equilibration_pbme_confined_cluster

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
   
   call do_shake(cluster,atomPairs,mdspecs)
   call get_distances_and_vectors(cluster,atomPairs)

   !call update_charges_in_complex_and_pairs(cluster,atomPairs)
   
   !get force
   call get_all_forces(atomPairs,force)

   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do

   call do_rattle(cluster,atomPairs,mdspecs)
end subroutine velocity_verlet_int_one_timestep

subroutine velocity_verlet_int_one_timestep_pbme(cluster,atomPairs,p,force,forceCCoM,mdspecs)
use stateevaluation
use forcecalculation
use quantumcalculations
use maproutines
implicit none
   real(8),dimension(1:3),intent(inout) :: forceCCoM
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: mdspecs
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j, nAtoms, nMap, nBasisFun
   real(8) :: dt,hdt
   real(8),dimension(1:3) :: CCoMvel
   real(8),dimension(:),allocatable :: allEigenVal
   real(8),dimension(:,:),allocatable :: ht, HMatrix, allEigenVec

   nAtoms = size(cluster)
   nMap = size(p%eigenvalues)
   nBasisFun = size(p%phi,1)
   
   allocate(ht(1:nMap,1:nMap))
   allocate(HMatrix(1:nBasisFun,1:nBasisFun))
   allocate(allEigenVec(1:nBasisFun,1:nBasisFun))
   allocate(allEigenVal(1:nBasisFun))

   dt = mdspecs%timeStep
   hdt = mdspecs%halfTimeStep
   
   !CCoMvel = forceToVelUnits*hdt*forceCCoM/(cluster(1)%mass+cluster(2)%mass)
   !cluster(1)%vel = cluster(1)%vel + CCoMVel*cluster(1)%mass/(cluster(1)%mass-cluster(2)%mass)
   !cluster(2)%vel = cluster(2)%vel - CCoMVel*cluster(2)%mass/(cluster(1)%mass-cluster(2)%mass)
   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
      cluster(j)%pos = cluster(j)%pos + dt*cluster(j)%vel
   end do
   
   do i = 1, nMap
      do j = 1, nMap
         p%pm(j) = p%pm(j) - (hdt/hbar)*p%h(j,i)*p%rm(i)
      end do
   end do 
   
   call do_shake(cluster,atomPairs,mdspecs)
   call get_distances_and_vectors(cluster,atomPairs)
   call get_H_grid_Atoms_pos_and_vec(p%gridHSolvent,atomPairs)
   
   !get new lambdas if requested
   !next line is necessary for h lambda lambda anyways
   call get_phi_Vsubsystem_phi_matrix(p%phi,atomPairs(1,2)%rij,p%phiVsphi)
   if (mdSpecs%updateLambdasOntheFly == 1) then
      HMatrix = p%phiKphi + p%phiVsphi
      call get_subsystem_lambdas(HMatrix,p%SMatrix,allEigenVec,allEigenVal)
      if (nMap > 2) then
         p%eigenvalues(1:nMap) = allEigenVal(1:nMap)
         p%lambda(1:nBasisFun,1:nMap) = allEigenVec(1:nBasisFun,1:nMap)
      else if (nMap == 2) then
         p%eigenvalues(1) = allEigenVal(1)
         p%eigenvalues(2) = allEigenVal(3)
         p%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,1)
         p%lambda(1:nBasisFun,2) = allEigenVec(1:nBasisFun,3)
      else if (nMap == 1) then
         if (mdSpecs%singleMap == 1) then
            p%eigenvalues(1) = allEigenVal(1)
            p%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,1)
         else
            p%eigenvalues(1) = allEigenVal(3)
            p%lambda(1:nBasisFun,1) = allEigenVec(1:nBasisFun,3)
         end if
      else
         stop 'error in number of quantum states classically mapped (check nMapStates)'
      end if
   end if
   
   !call necessary stuff to update h lambda lambda
   call get_mapFactor(p)
   call get_phi_inv_r_HS_phi_matrix(p)
   call get_lambda_h_lambda_matrix(cluster,atomPairs,p)
   call make_matrix_traceless(p%h,p%hTraceN,ht)
   p%h = ht

   do i = 1, nMap
      do j = 1, nMap
         p%rm(j) = p%rm(j) + (dt/hbar)*p%h(j,i)*p%pm(i)
      end do
   end do 
   
   do i = 1, nMap
      do j = 1, nMap
         p%pm(j) = p%pm(j) - (hdt/hbar)*p%h(j,i)*p%rm(i)
      end do
   end do 

   !call necessary stuff to update matrix elements to get force
   call get_phi_d_VBH_phi_matrix(p,atomPairs(1,2)%rij)
   !call get_phi_inv_r2_HS_phi_matrix(p)
   call get_phi_inv_r3_HS_phi_matrix(p)
   call get_phi_rc_inv_r3_HS_phi_matrix(atomPairs(1,2)%rij,p)
   call get_mapFactor(p)
   call get_all_forces_pbme(cluster,atomPairs,p,force,forceCCoM)

   !CCoMvel = forceToVelUnits*hdt*forceCCoM/(cluster(1)%mass+cluster(2)%mass)
   !cluster(1)%vel = cluster(1)%vel + CCoMVel*cluster(1)%mass/(cluster(1)%mass-cluster(2)%mass)
   !cluster(2)%vel = cluster(2)%vel - CCoMVel*cluster(2)%mass/(cluster(1)%mass-cluster(2)%mass)
   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do

   call do_rattle(cluster,atomPairs,mdspecs)
end subroutine velocity_verlet_int_one_timestep_pbme

subroutine velocity_verlet_int_one_timestep_pbme_confined_cluster(cluster,atomPairs,p,force,forceCCoM,mdspecs)
use stateevaluation
use forcecalculation
use quantumcalculations
use maproutines
implicit none
   real(8),dimension(1:3),intent(inout) :: forceCCoM
   type(Atom),dimension(:),intent(inout) :: cluster
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(inout) :: atomPairs
   type(MdData),intent(in) :: mdspecs
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j, nAtoms, nMap
   real(8) :: d,dt,hdt,clusterRadius, clusterMaxRadius
   real(8),dimension(1:3) :: CCoMvel,com
   real(8),dimension(:,:),allocatable :: ht

   nAtoms = size(cluster)
   nMap = size(p%eigenvalues)
   
   allocate(ht(1:nMap,1:nMap))

   dt = mdspecs%timeStep
   hdt = mdspecs%halfTimeStep
   
   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
      cluster(j)%pos = cluster(j)%pos + dt*cluster(j)%vel
   end do
   
   do i = 1, nMap
      do j = 1, nMap
         p%pm(j) = p%pm(j) - (hdt/hbar)*p%h(j,i)*p%rm(i)
      end do
   end do 
   
   call do_shake(cluster,atomPairs,mdspecs)
   call get_distances_and_vectors(cluster,atomPairs)
   call get_H_grid_Atoms_pos_and_vec(p%gridHSolvent,atomPairs)
   
   !call necessary stuff to update h lambda lambda
   call get_mapFactor(p)
   call get_phi_inv_r_HS_phi_matrix(p)
   call get_lambda_h_lambda_matrix(cluster,atomPairs,p)
   call make_matrix_traceless(p%h,p%hTraceN,ht)

   do i = 1, nMap
      do j = 1, nMap
         p%rm(j) = p%rm(j) + (dt/hbar)*p%h(j,i)*p%pm(i)
      end do
   end do 
   
   do i = 1, nMap
      do j = 1, nMap
         p%pm(j) = p%pm(j) - (hdt/hbar)*p%h(j,i)*p%rm(i)
      end do
   end do 

   !call necessary stuff to update matrix elements to get force
   call get_phi_d_VBH_phi_matrix(p,atomPairs(1,2)%rij)
   call get_phi_inv_r3_HS_phi_matrix(p)
   call get_phi_rc_inv_r3_HS_phi_matrix(atomPairs(1,2)%rij,p)
   call get_mapFactor(p)
   call get_all_forces_pbme(cluster,atomPairs,p,force,forceCCoM)

   !check if all atoms are in confinement otherwise add force to those not there
   clusterRadius = (nAtoms/(0.012d0))**(1d0/3d0)
   clusterMaxRadius = 2d0*clusterRadius
   call get_center_of_mass_vector(cluster,com)
   do i = 1, nAtoms
      d = sum((cluster(i)%pos - com)**2)
      if (d >= clusterMaxRadius) then
         force%inAtom(i)%total = force%inAtom(i)%total - (8.4d0*cluster(i)%pos - com)
      end if
   end do
   
   do j = 1, nAtoms
      cluster(j)%vel = cluster(j)%vel + forceToVelUnits*hdt*force%inAtom(j)%total/cluster(j)%mass
   end do

   call do_rattle(cluster,atomPairs,mdspecs)
end subroutine velocity_verlet_int_one_timestep_pbme_confined_cluster

subroutine do_shake(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md
   
   integer :: i,j,ai,aj
   real(8) :: esig1,amti,amtj,gamm,gammi,gammj
   real(8),dimension(1:3) :: tempBondVec
   real(8),dimension(:,:),allocatable :: originalBondVec,at_temp

   allocate(originalBondVec(1:md%nBondConstraints,1:3))
   allocate(at_temp(1:md%nBondConstraints*2,1:3))

   do i = 1, md%nBondConstraints
      ai = 3 + (i-1)*2
      aj = 4 + (i-1)*2
      originalBondVec(i,1:3) = pair(ai,aj)%vectorij

      at_temp(ai-2,1:3) = at(ai)%pos
      at_temp(aj-2,1:3) = at(aj)%pos
   end do
   
   do i = 1, md%nBondConstraints
      ai = 1 + (i-1)*2
      aj = 2 + (i-1)*2
      do j = 1, maxShakeCycles
         tempBondVec = at_temp(aj,1:3) - at_temp(ai,1:3)
         esig1 = sum(tempBondVec**2) - constrainedRS**2
         if (abs(esig1) < toleranceConstraints) exit

         amti = 1d0/at(ai+2)%mass
         amtj = 1d0/at(aj+2)%mass
         gamm = esig1/(2d0*(amti+amtj)*sum(tempBondVec*originalBondVec(i,1:3)))
         
         gammi = gamm*amti
         at_temp(ai,1:3) = at_temp(ai,1:3) + gammi*originalBondVec(i,1:3)
         at(ai+2)%vel = at(ai+2)%vel + gammi*originalBondVec(i,1:3)/md%timeStep
         gammj = gamm*amtj
         at_temp(aj,1:3) = at_temp(aj,1:3) - gammj*originalBondVec(i,1:3)
         at(aj+2)%vel = at(aj+2)%vel - gammj*originalBondVec(i,1:3)/md%timeStep
      end do
      if (j == maxShakeCycles) stop 'shake did not converge'
   end do

   do i = 1, md%nBondConstraints
      ai = 1 + (i-1)*2
      aj = 2 + (i-1)*2
      
      at(ai+2)%pos = at_temp(ai,1:3)
      at(aj+2)%pos = at_temp(aj,1:3)
   end do
end subroutine do_shake

subroutine do_shake_with_H(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md
   
   integer :: i,j,ai,aj
   real(8) :: esig1,amti,amtj,gamm,gammi,gammj
   real(8),dimension(1:3) :: tempBondVec
   real(8),dimension(:,:),allocatable :: originalBondVec,at_temp

   allocate(originalBondVec(1:md%nBondConstraints,1:3))
   allocate(at_temp(1:md%nBondConstraints*2,1:3))

   do i = 1, md%nBondConstraints
      ai = 4 + (i-1)*2
      aj = 5 + (i-1)*2
      originalBondVec(i,1:3) = pair(ai,aj)%vectorij

      at_temp(ai-3,1:3) = at(ai)%pos
      at_temp(aj-3,1:3) = at(aj)%pos
   end do
   
   do i = 1, md%nBondConstraints
      ai = 1 + (i-1)*2
      aj = 2 + (i-1)*2
      do j = 1, maxShakeCycles
         tempBondVec = at_temp(aj,1:3) - at_temp(ai,1:3)
         esig1 = sum(tempBondVec**2) - constrainedRS**2
         if (abs(esig1) < toleranceConstraints) exit

         amti = 1d0/at(ai+3)%mass
         amtj = 1d0/at(aj+3)%mass
         gamm = esig1/(2d0*(amti+amtj)*sum(tempBondVec*originalBondVec(i,1:3)))
         
         gammi = gamm*amti
         at_temp(ai,1:3) = at_temp(ai,1:3) + gammi*originalBondVec(i,1:3)
         at(ai+3)%vel = at(ai+3)%vel + gammi*originalBondVec(i,1:3)/md%timeStep
         gammj = gamm*amtj
         at_temp(aj,1:3) = at_temp(aj,1:3) - gammj*originalBondVec(i,1:3)
         at(aj+3)%vel = at(aj+3)%vel - gammj*originalBondVec(i,1:3)/md%timeStep
      end do
      if (j == maxShakeCycles) stop 'shake did not converge'
   end do

   do i = 1, md%nBondConstraints
      ai = 1 + (i-1)*2
      aj = 2 + (i-1)*2
      
      at(ai+3)%pos = at_temp(ai,1:3)
      at(aj+3)%pos = at_temp(aj,1:3)
   end do
end subroutine do_shake_with_H

subroutine do_rattle(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md

   integer :: i,j,ai,aj
   real(8) :: vvv,wwm,gmma
   
   do i = 1, md%nBondConstraints
      ai = 3 + (i-1)*2
      aj = 4 + (i-1)*2
      do j = 1, maxRattleCycles
         vvv = sum( pair(ai,aj)%vectorij*(at(aj)%vel - at(ai)%vel)  )
         if (abs(vvv) < toleranceConstraints) exit
         wwm = 1d0/at(ai)%mass + 1d0/at(aj)%mass
         gmma = vvv/(wwm*constrainedRS**2)!pair(ai,aj)%rij**2)
         
         at(ai)%vel = at(ai)%vel + gmma*pair(ai,aj)%vectorij/at(ai)%mass
         at(aj)%vel = at(aj)%vel - gmma*pair(ai,aj)%vectorij/at(aj)%mass
      end do
      if (j == maxShakeCycles) stop 'shake did not converge'
   end do
end subroutine do_rattle

subroutine do_rattle_with_H(at,pair,md)
implicit none
   type(Atom),dimension(:),intent(inout) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(MdData),intent(in) :: md

   integer :: i,j,ai,aj
   real(8) :: vvv,wwm,gmma
   
   do i = 1, md%nBondConstraints
      ai = 4 + (i-1)*2
      aj = 5 + (i-1)*2
      do j = 1, maxRattleCycles
         vvv = sum( pair(ai,aj)%vectorij*(at(aj)%vel - at(ai)%vel)  )
         if (abs(vvv) < toleranceConstraints) exit
         wwm = 1d0/at(ai)%mass + 1d0/at(aj)%mass
         gmma = vvv/(wwm*constrainedRS**2)!pair(ai,aj)%rij**2)
         
         at(ai)%vel = at(ai)%vel + gmma*pair(ai,aj)%vectorij/at(ai)%mass
         at(aj)%vel = at(aj)%vel - gmma*pair(ai,aj)%vectorij/at(aj)%mass
      end do
      if (j == maxShakeCycles) stop 'shake did not converge'
   end do
end subroutine do_rattle_with_H

end module dynamicsroutines
