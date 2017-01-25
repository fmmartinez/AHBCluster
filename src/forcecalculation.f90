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

subroutine get_all_forces_pbme(at,pairs,lambda,pAHp,pBHp,pqAp,pqBp,pir2p,mapFactor,force,forceAtComplexCoM)
!summary of cryptic names
!pir2p = <phi|1/r^2|phi> matrix (r is of H with all solvent atoms)
use quantumcalculations
use maproutines
implicit none
   
   type(Atom),dimension(:),intent(in) :: at
   type(Forces),intent(inout) :: force
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   real(8),dimension(1:3),intent(out) :: forceAtComplexCoM
   real(8),dimension(:,:),intent(in) :: lambda, pAHp, pBHp, pqAp, pqBp, mapFactor
   type(MatrixList),dimension(:,:),intent(in) :: pir2p
   
   integer :: i,j,n,nm
   type(VectorForMatrix),dimension(:),allocatable :: forceVecHAtomTemp
   real(8),dimension(:),allocatable :: forceHAtomTemp
   real(8),dimension(:,:),allocatable :: dh,dhx,dhy,dhz

   n = size(force%inAtom)
   nm = size(mapFactor,1)

   allocate(dh(1:nm,1:nm))
   allocate(dhx(1:nm,1:nm))
   allocate(dhy(1:nm,1:nm))
   allocate(dhz(1:nm,1:nm))
   allocate(forceHAtomTemp(1:n))
   allocate(forceVecHAtomTemp(1:n))

   do i = 1, n
      force%inAtom(i)%total = [0d0,0d0,0d0]
      forceVecHAtomTemp(i)%vecij = [0d0,0d0,0d0]
      forceHAtomTemp = 0d0
      do j = 1, n
         force%atomPair(i,j) = 0d0
      end do
   end do

   !complex forces   
   force%atomPair(1,2) = get_AB_force(pairs(1,2)%rij)
   force%atomPair(2,1) = force%atomPair(1,2)
      !these are forces due to the proton
   call get_lambda_d_VAH_lambda_matrix(lambda,pAHp,dh)
   forceHAtomTemp(1) = get_map_contribution(-dh,mapFactor)
   forceAtComplexCoM = -forceHAtomTemp(1)*pairs(2,1)%vectorij/pairs(2,1)%rij
   call get_lambda_d_VBH_lambda_matrix(lambda,pBHp,dh)
   forceHAtomTemp(2) = get_map_contribution(-dh,mapFactor)
   forceAtComplexCoM = forceAtComplexCoM + (-forceHAtomTemp(2)*pairs(1,2)%vectorij)/pairs(2,1)%rij
   
   !AB vs S
   do j = 3, n
      !A vs S
      i = 1
      force%atomPair(i,j) = get_lj_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%rij)
      call get_lambda_d_VASol_lambda_matrix(at(j),pairs(i,j),lambda,pqAp,dh)
      force%atomPair(i,j) = force%atomPair(i,j) + get_map_contribution(-dh,mapFactor)
      
      force%atomPair(j,i) = force%atomPair(i,j)
      
      !B vs S
      i = 2
      force%atomPair(i,j) = get_lj_force(pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%rij)
      call get_lambda_d_VBSol_lambda_matrix(at(j),pairs(i,j),lambda,pqBp,dh)
      force%atomPair(i,j) = force%atomPair(i,j) + get_map_contribution(-dh,mapFactor)
     
      force%atomPair(j,i) = force%atomPair(i,j)

      !H vs S
      !this is saved in a temporal variable
      call get_lambda_d_VHSol_lambda_matrix(at(j),lambda,pir2p(j,1:3),dhx,dhy,dhz)
      forceVecHAtomTemp(j)%vecij(1) = get_map_contribution(-dhx,mapFactor)
      forceVecHAtomTemp(j)%vecij(2) = get_map_contribution(-dhy,mapFactor)
      forceVecHAtomTemp(j)%vecij(3) = get_map_contribution(-dhz,mapFactor)
      forceHAtomTemp(j) = sum(forceVecHAtomTemp(j)%vecij**2)
      forceAtComplexCoM = forceAtComplexCoM + (-forceVecHAtomTemp(j)%vecij)
   end do
   
   !S vs S
   do i = 3, n-1, 2
      force%atomPair(i,i+1) = 0d0!get_SS_bond_force(pairs(i,i+1)%rij)
      force%atomPair(i+1,i) = force%atomPair(i,i+1)
      do j = i+2, n
         force%atomPair(i,j) = get_ljelec_force(&
                                 pairs(i,j)%ljEps,pairs(i,j)%ljSig,pairs(i,j)%qq,pairs(i,j)%rij)
         force%atomPair(i+1,j) = get_ljelec_force(&
                                 pairs(i+1,j)%ljEps,pairs(i+1,j)%ljSig,pairs(i+1,j)%qq,pairs(i+1,j)%rij)
         force%atomPair(j,i) = force%atomPair(i,j)
         force%atomPair(j,i+1) = force%atomPair(i+1,j)
      end do
   end do
   
   !add all forces due to interactions between atoms A, B and Solvent
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

   !add forces due to interactions of A/B/Solvent with H wavefunction
   !A vs H, unitary vector along A--B used
   force%inAtom(1)%total = force%inAtom(1)%total + &
                           forceHAtomTemp(1)*pairs(2,1)%vectorij/pairs(2,1)%rij
   !B vs H, unitary vector along B--A used
   force%inAtom(2)%total = force%inAtom(2)%total + &
                           forceHAtomTemp(2)*pairs(1,2)%vectorij/pairs(1,2)%rij
   !S vs H
   do i = 3, n
      force%inAtom(i)%total = force%inAtom(i)%total + forceVecHAtomTemp(i)%vecij
   end do
end subroutine get_all_forces_pbme

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

function get_ljelec_force(eps,sig,qq,r) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: eps,sig,qq,r

   f = 24d0*eps*(2d0*sig**12/r**13 - sig**6/r**7) + kCoulomb*qq/r**2

end function get_ljelec_force

function get_lj_force(eps,sig,r) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: eps,sig,r

   f = 24d0*eps*(2d0*sig**12/r**13 - sig**6/r**7)

end function get_lj_force

function get_elec_force(qq,r) result(f)
implicit none
   real(8) :: f
   real(8),intent(in) :: qq,r

   f = kCoulomb*qq/r**2

end function get_elec_force

end module forcecalculation
