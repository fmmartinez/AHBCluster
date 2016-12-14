module stateevaluation
use definitions

contains

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

subroutine update_charges_in_complex(atoms,pairs)
implicit none
   
   real(8) :: rah, fr
   type(Atom),dimension(:),intent(inout) :: atoms
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   
   rah = pairs(1,3)%rij
   fr = 0.5d0*(1d0 + (rah - 1.43d0)/sqrt((rah-1.43d0)**2+0.125d0**2))
   atoms(1)%charge = (1d0-fr)*(-0.5d0)+fr*(-1d0)
   atoms(2)%charge = fr*(0.5d0)

end subroutine update_charges_in_complex

end module stateevaluation
