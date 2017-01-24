module maproutines
use definitions

contains

function get_map_contribution(mat,mapFactor) result(y)
implicit none
   real(8),dimension(:,:),intent(in) :: mapFactor,mat
   
   integer :: i,j,nm
   real(8) :: y
   
   nm = size(mapFactor,1)

   y = 0d0
   do i = 1, nm
      do j = 1, nm
         y = y + mat(i,j)*mapFactor(i,j)
      end do
   end do
end function get_map_contribution

subroutine get_mapFactor(rm,pm,mapFactor)
implicit none
   real(8),dimension(:),intent(in) :: rm,pm
   real(8),dimension(:,:),intent(out) :: mapFactor
   
   integer :: i,j,nm

   nm = size(rm)

   mapFactor = 0d0
   do i = 1, nm
      do j = 1, nm
         mapFactor = 0.5d0*(rm(i)*rm(j) + pm(i)*pm(j))/hbar
      end do
   end do

end subroutine get_mapFactor

end module maproutines
