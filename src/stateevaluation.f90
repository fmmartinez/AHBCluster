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

end module stateevaluation
