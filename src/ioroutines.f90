module ioroutines
use definitions
implicit none

contains

subroutine read_md_input_file(n,nBasisFunCov,nBasisFunIon,nMapStates,md)
implicit none
   integer,intent(inout) :: n,nBasisFunCov,nBasisFunIon,nMapStates
   type(MdData),intent(inout) :: md

   integer :: unit1

   open(newunit=unit1,file='md.in')
      read(unit1,*)
      read(unit1,*) n, nBasisFunCov, nBasisFunIon, nMapStates, md%singleMap
      read(unit1,*)
      read(unit1,*) md%seed, md%nTrajectories
      read(unit1,*)
      read(unit1,*) md%timeStep
      read(unit1,*)
      read(unit1,*) md%eqSteps,md%maxEqTries,md%eqPhases,md%stepFreqEqSave
      read(unit1,*)
      read(unit1,*) md%initialEqTempInK,md%targetTempInK
      read(unit1,*)
      read(unit1,*) md%prodSteps
      read(unit1,*)
      read(unit1,*) md%stepFreqVelRescale,md%stepFreqCoMRemoval
      read(unit1,*)
      read(unit1,*) md%stepFreqOutTrajectory, md%stepFreqOutLog
      read(unit1,*)
      read(unit1,*) md%updateLambdasOntheFly
   close(unit1)
   
   md%halfTimeStep = md%timeStep/2d0

   md%eqPhaseSteps = md%eqSteps/md%eqPhases

   md%nBondConstraints = (n-2)/2
end subroutine read_md_input_file

subroutine read_force_field_file(at)
implicit none
   type(atom),dimension(:),intent(inout) :: at

   integer :: i,unit1,n

   n = size(at)
   
   open(newunit=unit1,file='field.in')
   read(unit1,*) 
   do i = 1, n
      read(unit1,'(a4,4f14.8)') at(i)%symbol, at(i)%mass, at(i)%charge, &
                                 at(i)%ljSigma, at(i)%ljEpsilon
   end do
   close(unit1)

end subroutine read_force_field_file

subroutine initialize_force_field_no_H(at)
   type(atom),dimension(:),intent(inout) :: at
   
   integer :: i,n

   n = size(at)
   if (mod(n,2) /= 0) stop 'number of atoms must be even if no H is present'

   at(1)%symbol = 'O'
   at(1)%mass = 93d0
   at(1)%charge = 0.0d0 !placeholder, it is replaced latter
   at(1)%ljSigma = 3.5d0
   at(1)%ljEpsilon = 0.3974d0

   at(2)%symbol = 'N'
   at(2)%mass = 59d0
   at(2)%charge = 0d0 !placeholder, it is replaced latter
   at(2)%ljSigma = 3.5d0
   at(2)%ljEpsilon = 0.3974d0

   do i = 3, n, 2
      at(i)%symbol = 'C'
      at(i)%mass = 15d0
      at(i)%charge = 0.25d0
      at(i)%ljSigma = 3.774d0
      at(i)%ljEpsilon = 0.238d0
      
      at(i+1)%symbol = 'Cl'
      at(i+1)%mass = 35.5d0
      at(i+1)%charge = -0.25d0
      at(i+1)%ljSigma = 3.481d0
      at(i+1)%ljEpsilon = 0.4150d0
   end do
end subroutine initialize_force_field_no_H

subroutine initialize_force_field_explicit_H(at)
   type(atom),dimension(:),intent(inout) :: at
   
   integer :: i,n

   n = size(at)

   at(1)%symbol = 'O'
   at(1)%mass = 93d0
   at(1)%charge = -0.5d0 !placeholder, it is replaced latter
   at(1)%ljSigma = 3.5d0
   at(1)%ljEpsilon = 0.3974d0

   at(2)%symbol = 'N'
   at(2)%mass = 59d0
   at(2)%charge = 0d0 !placeholder, it is replaced latter
   at(2)%ljSigma = 3.5d0
   at(2)%ljEpsilon = 0.3974d0

   at(3)%symbol = 'H'
   at(3)%mass = 1d0
   at(3)%charge = 0.5d0
   at(3)%ljSigma = 0d0
   at(3)%ljEpsilon = 0d0

   do i = 4, n, 2
      at(i)%symbol = 'C'
      at(i)%mass = 15d0
      at(i)%charge = 0.25d0
      at(i)%ljSigma = 3.774d0
      at(i)%ljEpsilon = 0.238d0
      
      at(i+1)%symbol = 'Cl'
      at(i+1)%mass = 35.5d0
      at(i+1)%charge = -0.25d0
      at(i+1)%ljSigma = 3.481d0
      at(i+1)%ljEpsilon = 0.4150d0
   end do
end subroutine initialize_force_field_explicit_H

subroutine read_config_in_XYZ_file(at)
implicit none
   type(atom),dimension(:),intent(inout) :: at

   integer :: i,unit1,n

   n = size(at)
        
   open(newunit=unit1,file='config_in.xyz')
   read(unit1,*)
   read(unit1,*)
   do i = 1, n
      read(unit1,'(a4,3f14.8)') at(i)%symbol, at(i)%pos(1:3)
   end do
   close(unit1)

end subroutine read_config_in_XYZ_file

subroutine write_generated_initial_positions_XYZ(at)
implicit none
   type(atom),dimension(:),intent(in) :: at
   
   integer :: unit1,nAtoms,j

   nAtoms = size(at)

   open(newunit=unit1,file='ini_generated.xyz')
   write(unit1,*) nAtoms
   write(unit1,*) 'symbol - positions x y z -- intial config'
   do j = 1, nAtoms
      write(unit1,'(a4,3f14.8)') at(j)%symbol, at(j)%pos(1:3)
   end do
   close(unit1)
end subroutine write_generated_initial_positions_XYZ

subroutine write_xyz_trajectory(at,step,unitNumber)
implicit none
   integer,intent(in) :: step,unitNumber
   type(atom),dimension(:),intent(in) :: at
   
   integer :: nAtoms,j

   nAtoms = size(at)
   
   write(unitNumber,*) nAtoms
   write(unitNumber,*) 'symbol - positions x y z --', step
   do j = 1, nAtoms
      write(unitNumber,'(a4,3f14.8)') at(j)%symbol, at(j)%pos(1:3)
   end do
end subroutine write_xyz_trajectory

subroutine write_equilibration_trajectory(at,i)
implicit none
   integer,intent(in) :: i
   type(atom),dimension(:),intent(in) :: at
   
   integer :: unit1,nAtoms,j

   nAtoms = size(at)
   
   open(newunit=unit1,file='trajectory-eq.xyz',position='append')
   write(unit1,*) nAtoms
   write(unit1,*) 'symbol - positions x y z --',i
   do j = 1, nAtoms
      write(unit1,'(a4,3f14.8)') at(j)%symbol, at(j)%pos(1:3)
   end do
   close(unit1)
end subroutine write_equilibration_trajectory

subroutine write_production_trajectory(at,i)
implicit none
   integer,intent(in) :: i
   type(atom),dimension(:),intent(in) :: at
   
   integer :: unit1,nAtoms,j

   nAtoms = size(at)
   
   open(newunit=unit1,file='trajectory-pd.xyz',position='append')
   write(unit1,*) nAtoms
   write(unit1,*) 'symbol - positions x y z --',i
   do j = 1, nAtoms
      write(unit1,'(a4,3f14.8)') at(j)%symbol, at(j)%pos(1:3)
   end do
   close(unit1)
end subroutine write_production_trajectory

end module ioroutines
