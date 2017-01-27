module stateevaluation
use definitions

contains

subroutine generate_positions(atoms)
implicit none
   
   integer :: i,j,n,nmol,l,counter
   
   real(8),dimension(1:12,1:3) :: v
   real(8),dimension(:,:),allocatable :: cen

   type(Atom),dimension(:),intent(inout) :: atoms
   
   !fcc generating vectors
   v(1,1:3) = [0,1,1]
   v(2,1:3) = [1,0,1]
   v(3,1:3) = [1,1,0]

   v(4,1:3) = [0,-1,1]
   v(5,1:3) = [-1,0,1]
   v(6,1:3) = [-1,1,0]
   
   v(7,1:3) = [0,1,-1]
   v(8,1:3) = [1,0,-1]
   v(9,1:3) = [1,-1,0]
   
   v(10,1:3) = [0,-1,-1]
   v(11,1:3) = [-1,0,-1]
   v(12,1:3) = [-1,-1,0]
   v = v/sqrt(2d0)


   n = size(atoms)
   nmol = (n - 2)/2 + 1
   l = nmol/12
   
   !add geometrical centers
   allocate(cen(1:nmol,1:3))
   counter = 0d0
   do i = 1, l+1
      do j = 1, 12
         counter = counter + 1
         if (counter == nmol) exit
         cen(j+(i*12-12),1:3) = v(j,1:3)*5.3d0*i
      end do
   end do
   
   !place complex
   atoms(1)%pos = cen(1,1:3) + [-1.35d0,0.0d0,0.0d0]
   atoms(2)%pos = cen(1,1:3) + [ 1.35d0,0.0d0,0.0d0]

   !place solvent
   do i = 2, nmol
      atoms(i*2-1)%pos = cen(i,1:3) + [-0.89d0,0d0,0d0]
      atoms(i*2)%pos = cen(i,1:3)  + [0.89d0,0d0,0d0]
   end do
end subroutine generate_positions

subroutine generate_positions_with_H(atoms)
implicit none
   
   integer :: i,j,n,nmol,l,counter
   
   real(8),dimension(1:12,1:3) :: v
   real(8),dimension(:,:),allocatable :: cen

   type(Atom),dimension(:),intent(inout) :: atoms
   
   !fcc generating vectors
   v(1,1:3) = [0,1,1]
   v(2,1:3) = [1,0,1]
   v(3,1:3) = [1,1,0]

   v(4,1:3) = [0,-1,1]
   v(5,1:3) = [-1,0,1]
   v(6,1:3) = [-1,1,0]
   
   v(7,1:3) = [0,1,-1]
   v(8,1:3) = [1,0,-1]
   v(9,1:3) = [1,-1,0]
   
   v(10,1:3) = [0,-1,-1]
   v(11,1:3) = [-1,0,-1]
   v(12,1:3) = [-1,-1,0]
   v = v/sqrt(2d0)


   n = size(atoms)
   nmol = (n - 3)/2 + 1
   l = nmol/12
   
   !add geometrical centers
   allocate(cen(1:nmol,1:3))
   counter = 0d0
   do i = 1, l+1
      do j = 1, 12
         cen(j+(i*12-12),1:3) = v(j,1:3)*5.3d0*i
         counter = counter + 1
         if (counter == nmol) exit
      end do
   end do
   
   !place complex
   atoms(1)%pos = cen(1,1:3) + [-1.35d0,0.0d0,0.0d0]
   atoms(2)%pos = cen(1,1:3) + [ 1.35d0,0.0d0,0.0d0]
   atoms(3)%pos = cen(1,1:3) + [-0.35d0,0.0d0,0.0d0]

   !place solvent
   do i = 2, nmol
      atoms(i*2)%pos = cen(i,1:3) + [-0.89d0,0d0,0d0]
      atoms(i*2+1)%pos = cen(i,1:3)  + [0.89d0,0d0,0d0]
   end do
   
   !m = n
   !!add complex so half lies on 0,0,0
   !atoms(1)%pos = [-1.35d0,0.0d0,0.0d0]
   !atoms(2)%pos = [ 1.35d0,0.0d0,0.0d0]
   !atoms(3)%pos = [-0.35d0,0.0d0,0.0d0]

   !!add centers so 

   !!do i = 4, n, 2
   !i = 4
   !counter = 0
   !do 
   !   vec = [ran(n)-0.5d0,ran(n)-0.5d0,ran(n)-0.5d0]
   !   dvec = sqrt(sum(vec**2))
   !   uvec = vec/dvec 
   !   atoms(i)%pos = uvec*(6d0*ran(n)) + [1.35d0,0d0,0d0]

   !   vec = [ran(n)-0.5d0,ran(n)-0.5d0,ran(n)-0.5d0]
   !   dvec = sqrt(sum(vec**2))
   !   uvec = vec/dvec
   !   atoms(i+1)%pos = 1.781d0*uvec + atoms(i)%pos

   !   do j = 1, i-1
   !      r1 = sqrt(sum( (atoms(i)%pos - atoms(j)%pos)**2 ))
   !      r2 = sqrt(sum( (atoms(i+1)%pos - atoms(j)%pos)**2 ))
   !      if (r1<=3.5 .or. r2<=3.5) counter = counter + 1
   !   end do

   !   if (counter > 0) then
   !      counter = 0
   !   else
   !      counter = 0
   !      i = i + 2
   !   end if

   !   if (i > m) exit
   !end do
end subroutine generate_positions_with_H

subroutine get_force_field_pair_parameters(atoms,pairs)
implicit none
   
   integer :: i,j,n

   type(Atom),dimension(:),intent(in) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)
   !initialize
   do i = 1, n
      do j = 1, n
         pairs(i,j)%ljEps = 0d0
         pairs(i,j)%ljSig = 0d0
         pairs(i,j)%qq = 0d0
      end do
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%qq = atoms(i)%charge*atoms(j)%charge
         pairs(j,i)%qq = pairs(i,j)%qq
      end do
   end do

   do j = 3, n
      i = 1
      pairs(i,j)%ljEps = atoms(i)%ljEpsilon
      pairs(j,i)%ljEps = pairs(i,j)%ljEps
      pairs(i,j)%ljSig = atoms(i)%ljSigma
      pairs(j,i)%ljSig = pairs(i,j)%ljSig
      
      i = 2
      pairs(i,j)%ljEps = atoms(i)%ljEpsilon
      pairs(j,i)%ljEps = pairs(i,j)%ljEps
      pairs(i,j)%ljSig = atoms(i)%ljSigma
      pairs(j,i)%ljSig = pairs(i,j)%ljSig
   end do

   do i = 3, n-1
      do j = i+1, n
         pairs(i,j)%ljEps = sqrt(atoms(i)%ljEpsilon*atoms(j)%ljEpsilon)
         pairs(j,i)%ljEps = pairs(i,j)%ljEps
         pairs(i,j)%ljSig = sqrt(atoms(i)%ljSigma*atoms(j)%ljSigma)
         pairs(j,i)%ljSig = pairs(i,j)%ljSig
      end do
   end do
end subroutine get_force_field_pair_parameters

subroutine get_force_field_pair_parameters_with_H(atoms,pairs)
implicit none
   
   integer :: i,j,n

   type(Atom),dimension(:),intent(in) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)
   !initialize
   do i = 1, n
      do j = 1, n
         pairs(i,j)%ljEps = 0d0
         pairs(i,j)%ljSig = 0d0
         pairs(i,j)%qq = 0d0
      end do
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%qq = atoms(i)%charge*atoms(j)%charge
         pairs(j,i)%qq = pairs(i,j)%qq
      end do
   end do

   do j = 4, n
      i = 1
      !pairs(i,j)%ljEps = sqrt(atoms(i)%ljEpsilon*atoms(j)%ljEpsilon)
      pairs(i,j)%ljEps = atoms(i)%ljEpsilon
      pairs(j,i)%ljEps = pairs(i,j)%ljEps
      !pairs(i,j)%ljSig = (atoms(i)%ljSigma + atoms(j)%ljSigma)/2
      pairs(i,j)%ljSig = atoms(i)%ljSigma
      pairs(j,i)%ljSig = pairs(i,j)%ljSig
      i = 2
      !pairs(i,j)%ljEps = sqrt(atoms(i)%ljEpsilon*atoms(j)%ljEpsilon)
      pairs(i,j)%ljEps = atoms(i)%ljEpsilon
      pairs(j,i)%ljEps = pairs(i,j)%ljEps
      !pairs(i,j)%ljSig = (atoms(i)%ljSigma + atoms(j)%ljSigma)/2
      pairs(i,j)%ljSig = atoms(i)%ljSigma
      pairs(j,i)%ljSig = pairs(i,j)%ljSig
      i = 3
      pairs(i,j)%ljEps = 0d0
      pairs(j,i)%ljEps = pairs(i,j)%ljEps
      pairs(i,j)%ljSig = 0d0
      pairs(j,i)%ljSig = pairs(i,j)%ljSig
   end do

   do i = 4, n-1
      do j = i+1, n
         pairs(i,j)%ljEps = sqrt(atoms(i)%ljEpsilon*atoms(j)%ljEpsilon)
         pairs(j,i)%ljEps = pairs(i,j)%ljEps
         !pairs(i,j)%ljSig = (atoms(i)%ljSigma + atoms(j)%ljSigma)/2
         pairs(i,j)%ljSig = sqrt(atoms(i)%ljSigma*atoms(j)%ljSigma)
         pairs(j,i)%ljSig = pairs(i,j)%ljSig
      end do
   end do
end subroutine get_force_field_pair_parameters_with_H

subroutine get_distances_and_vectors(atoms,pairs)
implicit none
   
   integer :: i,j,n

   type(Atom),dimension(:),intent(in) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)

   do i = 1, n
      pairs(i,i)%vectorij = [0d0,0d0,0d0]
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%vectorij = atoms(j)%pos - atoms(i)%pos
         pairs(j,i)%vectorij = -pairs(i,j)%vectorij
      end do
   end do
   
   do i = 1, n
      pairs(i,i)%rij = 0d0
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%rij = sqrt(sum(pairs(i,j)%vectorij**2))
         pairs(j,i)%rij = pairs(i,j)%rij
      end do
   end do

end subroutine get_distances_and_vectors

subroutine get_H_grid_Atoms_pos_and_vec(HS,pair)
!vectors have direction from atom to point in grid
implicit none
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(EvalOnGridHData),dimension(:),intent(out) :: HS
   
   integer :: i,j,n
   real(8) :: q
   
   n = size(pair,1)

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth

      HS(1)%gridPoint(i)%rij = q
      HS(1)%gridPoint(i)%vectorij = q*pair(1,2)%vectorij/pair(1,2)%rij
      
      HS(2)%gridPoint(i)%rij = pair(1,2)%rij - q
      HS(2)%gridPoint(i)%vectorij = HS(2)%gridPoint(i)%rij*pair(2,1)%vectorij/pair(1,2)%rij
      
      do j = 3, n
         HS(j)%gridpoint(i)%vectorij = pair(j,1)%vectorij + HS(1)%gridPoint(i)%vectorij
         HS(j)%gridpoint(i)%rij = sqrt(sum(HS(j)%gridpoint(i)%vectorij**2))
      end do
   end do
end subroutine get_H_grid_Atoms_pos_and_vec

subroutine update_charges_in_complex_and_pairs(atoms,pairs)
implicit none
   
   integer :: i,j,n

   real(8) :: rah, fr
   
   type(Atom),dimension(:),intent(inout) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)

   rah = pairs(1,3)%rij
   fr = 0.5d0*(1d0 + (rah - 1.43d0)/sqrt((rah-1.43d0)**2+0.125d0**2))
   
   atoms(1)%charge = (1d0-fr)*(-0.5d0)+fr*(-1d0)
   atoms(2)%charge = fr*(0.5d0)
   
   do i = 1, 2
      do j = 4, n
         pairs(i,j)%qq = atoms(i)%charge*atoms(j)%charge
         pairs(j,i)%qq = pairs(i,j)%qq
      end do
   end do

end subroutine update_charges_in_complex_and_pairs

subroutine get_center_of_mass_vector(at,com)
implicit none
   
   type(Atom),dimension(:),intent(in) :: at
   real(8),dimension(1:3),intent(out) :: com

   integer :: i,n
   real(8) :: totalMass

   n = size(at)
               
   totalMass = sum(at(1:n)%mass)

   com = [0d0,0d0,0d0]
   do i = 1, n
      com = com + at(i)%mass*at(i)%pos/totalMass
   end do
end subroutine get_center_of_mass_vector

function get_distance_solvent_CoM_complex_CoM(at) result(d)
implicit none

   type(Atom),dimension(:),intent(in) :: at
   
   integer :: n

   real(8) :: d
   real(8),dimension(1:3) :: complexCoM,solventCoM,dvec
   
   n = size(at)

   call get_center_of_mass_vector(at(1:2),complexCoM)
   call get_center_of_mass_vector(at(3:n),solventCoM)

   dvec = complexCom - solventCoM

   d = sqrt(sum(dvec**2))

end function get_distance_solvent_CoM_complex_CoM

function get_distance_solvent_CoM_complex_CoM_with_H(at) result(d)
implicit none

   type(Atom),dimension(:),intent(in) :: at
   
   integer :: n

   real(8) :: d
   real(8),dimension(1:3) :: complexCoM,solventCoM,dvec
   
   n = size(at)

   call get_center_of_mass_vector(at(1:3),complexCoM)
   call get_center_of_mass_vector(at(4:n),solventCoM)

   dvec = complexCom - solventCoM

   d = sqrt(sum(dvec**2))

end function get_distance_solvent_CoM_complex_CoM_with_H

subroutine generate_velocities(at,stream,tempInK)
use mkl_vsl_type
use mkl_vsl
implicit none

   real(8),intent(in) :: tempInK
   type(vsl_stream_state),intent(in) :: stream
   type(Atom),dimension(:),intent(inout) :: at
         
   integer :: i,errcode,nAtoms
   real(8),dimension(:),allocatable :: velSigma
   
   nAtoms = size(at)

   allocate(velSigma(1:nAtoms))

   velSigma = sqrt((tempInK/at%mass)*kBoltzmann)
                     
   do i = 1, nAtoms
      errcode = vdrnggaussian(method,stream,1,at(i)%vel(1),0d0,velSigma(i))
      errcode = vdrnggaussian(method,stream,1,at(i)%vel(2),0d0,velSigma(i))
      errcode = vdrnggaussian(method,stream,1,at(i)%vel(3),0d0,velSigma(i))
   end do

   deallocate(velSigma)
end subroutine generate_velocities

subroutine do_velocity_rescale(at,tempInK,nc)
implicit none
   integer,intent(in) :: nc
   real(8),intent(in) :: tempInK
   type(Atom),dimension(:),intent(inout) :: at

   integer :: i,nAtoms
   real(8) :: kineticEnergyInsta,kineticEnergyDesired,scaleFactor
   
   nAtoms = size(at)

   kineticEnergyInsta = get_kinetic_energy(at)/KToKcalMol

   kineticEnergyDesired = 0.5d0*(3d0*nAtoms-nc)*kBoltzmann*tempInK

   scaleFactor = sqrt(kineticEnergyDesired/kineticEnergyInsta)
   do i = 1, nAtoms
      at(i)%vel = scaleFactor*at(i)%vel
   end do
end subroutine do_velocity_rescale

function get_kinetic_energy(at) result (k)
implicit none
type(Atom),dimension(:),intent(in) :: at

   integer :: i,n
   real(8) :: k

   n = size(at)

   k = 0d0
   do i = 1, n
      k = k + sum(at(i)%mass*at(i)%vel**2)
   end do
   k = 0.5d0*KtoKcalMol*k
end function get_kinetic_energy

subroutine remove_CoM_movement(atoms)
implicit none
   
   integer :: i,n
   real(8) :: totalMass
   real(8),dimension(1:3) :: CoMVel

   type(Atom),dimension(:),intent(inout) :: atoms
   
   n = size(atoms)

   totalMass = sum(atoms(1:n)%mass)

   CoMVel = [0d0,0d0,0d0]
   do i = 1, n
      CoMVel = CoMVel + atoms(i)%mass*atoms(i)%vel/totalMass
   end do

   do i = 1, n
      atoms(i)%vel = atoms(i)%vel - CoMVel/n
   end do

end subroutine remove_CoM_movement

function get_total_momentum_magnitude(at) result (p)
   implicit none
      
   type(Atom),dimension(:),intent(in) :: at

   integer :: i,nAtoms
   real(8) :: p
   real(8),dimension(1:3) :: pVec
   
   nAtoms = size(at)

   pVec= 0d0
   do i = 1, nAtoms
      pVec = pVec + at(i)%mass*at(i)%vel
   end do

   p = sqrt(sum(pVec**2))
   end function get_total_momentum_magnitude

function get_insta_temperature(k,nAtoms,nc) result(t)
implicit none
   integer,intent(in) :: nAtoms,nc
   real(8),intent(in) :: k

   real(8) :: t

   t = (2d0/(3d0*nAtoms-nc))*k/ktoKcalMol
   t = t/kBoltzmann
end function get_insta_temperature

function get_solvent_polarization(at,pairs) result(de)
implicit none
   real(8),parameter :: s1 = 0d0, s2 = -0.56d0

   type(Atom),dimension(:),intent(in) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   
   integer :: i,nAtoms
   real(8) :: de,dSol1,dSol2
   real(8),dimension(1:3) :: com

   nAtoms = size(at)
   
   call get_center_of_mass_vector(at(1:2),com)
   de = 0d0
   do i = 3, nAtoms
      dSol1 = sqrt(sum(at(i)%pos - com)**2)
      dSol2 = sqrt(sum(at(i)%pos - (com - s2*pairs(1,2)%vectorij/pairs(1,2)%rij))**2)
      de = de + at(i)%charge*(1d0/dSol1 - 1d0/dSol2)
   end do
end function get_solvent_polarization

end module stateevaluation
