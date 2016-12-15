module stateevaluation
use definitions

contains

subroutine get_force_field_pair_parameters(atoms,pairs)
implicit none
   
   integer :: i,j,n

   type(Atom),dimension(:),intent(in) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%qq = atoms(i)%charge*atoms(j)%charge
         pairs(j,i)%qq = pairs(i,j)%qq

         pairs(i,j)%ljEps = sqrt(atoms(i)%ljEpsilon*atoms(j)%ljEpsilon)
         pairs(j,i)%ljEps = pairs(i,j)%ljEps

         pairs(i,j)%ljSig = (atoms(i)%ljSigma + atoms(j)%ljSigma)/2d0
         pairs(j,i)%ljSig = pairs(i,j)%ljSig
      end do
   end do
end subroutine get_force_field_pair_parameters

subroutine get_distances_and_vectors(atoms,pairs)
implicit none
   
   integer :: i,j,n

   type(Atom),dimension(:),intent(in) :: atoms
   type(AtomPairData),dimension(:,:),intent(inout) :: pairs
   
   n = size(atoms)

   do i = 1, n
      pairs(i,i)%rij = 0d0
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%rij = sqrt(sum((atoms(i)%pos - atoms(j)%pos)**2))
         pairs(j,i)%rij = pairs(i,j)%rij
      end do
   end do

   do i = 1, n
      pairs(i,i)%vectorij = [0d0,0d0,0d0]
   end do

   do i = 1, n-1
      do j = i+1, n
         pairs(i,j)%vectorij = atoms(j)%pos - atoms(i)%pos
         pairs(j,i)%vectorij = -pairs(i,j)%vectorij
      end do
   end do
end subroutine get_distances_and_vectors

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

end module stateevaluation
