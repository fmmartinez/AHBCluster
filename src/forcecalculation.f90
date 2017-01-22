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

   !complex forces   
   force%atomPair(1,2) = get_AB_force(pairs(1,2)%rij) + get_BH_force(pairs(1,2)%rij-1d0)
   force%atomPair(2,1) = force%atomPair(1,2)

   !AHB vs S
   do j = 3, n
      i = 1
      force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
      force%atomPair(j,i) = force%atomPair(i,j)
      
      i = 2
      force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
      force%atomPair(j,i) = force%atomPair(i,j)
   end do
   
   !S vs S
   do i = 3, n-1, 2
      force%atomPair(i,i+1) = 0d0!get_SS_bond_force(pairs(i,i+1)%rij)
      force%atomPair(i+1,i) = force%atomPair(i,i+1)
      do j = i+2, n
         force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
         force%atomPair(i+1,j) = get_ljelec_force(pairs(i+1,j)%ljEps,pairs(i+1,j)%ljSig,pairs(i+1,j)%qq,pairs(i+1,j)%rij)
         force%atomPair(j,i) = force%atomPair(i,j)
         force%atomPair(j,i+1) = force%atomPair(i+1,j)
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

subroutine get_all_forces_with_H(pairs,force)
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

   !complex forces   
   force%atomPair(1,2) = get_AB_force(pairs(1,2)%rij)
   force%atomPair(2,1) = force%atomPair(1,2)

   force%atomPair(1,3) = get_AH_force(pairs(1,3)%rij)
   force%atomPair(3,1) = force%atomPair(1,3)
   
   force%atomPair(2,3) = get_BH_force(pairs(2,3)%rij)
   force%atomPair(3,2) = force%atomPair(2,3)
   
   !AHB vs S
   do j = 4, n
      i = 1
      force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
      force%atomPair(j,i) = force%atomPair(i,j)
      
      i = 2
      force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
      force%atomPair(j,i) = force%atomPair(i,j)
      
      i = 3
      force%atomPair(i,j) = get_HS_force(pairs(i,j)%qq,pairs(i,j)%rij)
      force%atomPair(j,i) = force%atomPair(i,j)
   end do
   
   !S vs S
   do i = 4, n-1, 2
      force%atomPair(i,i+1) = 0d0!get_SS_bond_force(pairs(i,i+1)%rij)
      force%atomPair(i+1,i) = force%atomPair(i,i+1)
      do j = i+2, n
         force%atomPair(i,j) = get_ljelec_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
         force%atomPair(i+1,j) = get_ljelec_force(pairs(i+1,j)%ljEps,pairs(i+1,j)%ljSig,pairs(i+1,j)%qq,pairs(i+1,j)%rij)
         force%atomPair(j,i) = force%atomPair(i,j)
         force%atomPair(j,i+1) = force%atomPair(i+1,j)
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

end subroutine get_all_forces_with_H

function get_SS_bond_force(rij) result (f)
implicit none
   real(8),parameter :: kBond = 150d0   

   real(8) :: f
   real(8),intent(in) :: rij
   
   f = -2d0*kBond*(rij - 1.781d0)
end function get_SS_bond_force

function get_HS_force(qq,r) result (f)
implicit none
   real(8) :: f
   real(8),intent(in) :: qq,r

   f = kCoulomb*qq/r**2

end function get_HS_force

function get_ljelec_force(eps,sig,qq,r) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: eps,sig,qq,r

   f = 24d0*eps*(2d0*sig**12/r**13 - sig**6/r**7) + kCoulomb*qq/r**2

end function get_ljelec_force

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
