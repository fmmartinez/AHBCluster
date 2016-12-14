module forcecalculation
use definitions
implicit none

contains

subroutine get_all_forces(pairs,force)
   implicit none
   
   integer :: i,j,n

   type(Forces) :: force
   type(AtomPairData),dimension(:,:) :: pairs
   
   n = size(force%inAtom)

   do i = 1, n
      force%inAtom(i)%total = [0d0,0d0,0d0]
      do j = 1, n
         force%atomPair(i,j) = 0d0
      end do
   end do

   do i = 1, n-1
      do j = i+1, n
         if (i==1.and.j==2) force%atomPair(i,j) = get_AB_force(pairs(i,j)%rij)
         if (i==1.and.j==3) force%atomPair(i,j) = get_AH_force(pairs(i,j)%rij)
         !if (i==1.and.j>3) force%atomPair(i,j) = get_AS_force()
         if (i==2.and.j==3) force%atomPair(i,j) = get_BH_force(pairs(i,j)%rij)
         !if (i==2.and.j>3) force%atomPair(i,j) = get_BS_force()
         !if (i==3.and.j>3) force%atomPair(i,j) = get_SS_force()
         force%atomPair(j,i) = force%atomPair(i,j)
      end do 
   end do
   
   do i = 1, n
      do j = 1, i-1
         force%inAtom(i)%total = force%inAtom(i)%total + &
                                 force%atomPair(i,j)*pairs(j,i)%vectorij/pairs(j,i)%rij
      end do
      do j = i+1, n
         force%inAtom(i)%total = force%inAtom(i)%total + &
                                 force%atomPair(i,j)*pairs(j,i)%vectorij/pairs(j,i)%rij
      end do
   end do

end subroutine get_all_forces

function get_BH_force(rbh) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: rbh

   f = c*d*exp(-nb*(rbh-db)**2/(2d0*rbh))*(nb*(rbh-db)/rbh)*((rbh-db)/(2d0*rbh)-1d0)

end function get_BH_force

function get_AH_force(rah) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: rah

   f = d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)

end function get_AH_force

function get_AB_force(rab) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: rab

   f = (a*b)*exp(-a*rab)

end function get_AB_force

end module forcecalculation
