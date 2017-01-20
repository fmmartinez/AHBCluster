module ioroutines
use definitions
implicit none

contains

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

end module ioroutines
