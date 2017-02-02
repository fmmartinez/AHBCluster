module energycalculation
use definitions
implicit none

contains

function get_complex_energy(rah,rab) result(e)
implicit none
   real(8),intent(in) :: rab,rah

   real(8) :: e,rbh

   rbh = rab - rah
   e = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rbh-db)**2/(2d0*rbh)))
end function get_complex_energy

function get_complex_energy_repulsion_part(rab) result(e)
implicit none
   real(8),intent(in) :: rab

   real(8) :: e

   e = b*exp(-a*rab)

end function get_complex_energy_repulsion_part

function get_complex_energy_attraction_part(rah,rab) result(e)
implicit none
   real(8),intent(in) :: rab,rah

   real(8) :: e,rbh

   rbh = rab - rah
   e = d*(1-exp(-(na*(rah-da)**2)/(2d0*rah))) + c*d*(1-exp(-nb*(rbh-db)**2/(2d0*rbh)))

end function get_complex_energy_attraction_part

function get_complex_energy_with_H(pairs) result(e)
implicit none
   real(8) :: e
   real(8) :: rab,rah,rbh
   type(AtomPairData),dimension(:,:),intent(in) :: pairs

   rab = pairs(1,2)%rij
   rah = pairs(1,3)%rij
   rbh = pairs(2,3)%rij
   e = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rbh-db)**2/(2d0*rbh)))
end function get_complex_energy_with_H

subroutine get_complexsolvent_energy(pairs,elj,eel,e)
implicit none

   integer :: i,j,n
   real(8),intent(out) :: elj,eel,e
   real(8) :: r,sig,eps,qq
   type(AtomPairData),dimension(:,:),intent(in) :: pairs

   n = size(pairs,1)
   
   elj = 0d0   
   do i = 1, 2
      do j = 3, n
         r = pairs(i,j)%rij
         sig = pairs(i,j)%ljSig
         eps = pairs(i,j)%ljEps
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
      end do
   end do
   
   eel = 0d0   
   do i = 1, 2 
      do j = 3, n
         r = pairs(i,j)%rij
         qq = pairs(i,j)%qq
         eel = eel + kCoulomb*qq/r
      end do
   end do
   
   e = elj + eel   
end subroutine get_complexsolvent_energy

function get_complexsolvent_lj_energy(pairs) result(elj)
implicit none
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   
   integer :: i,j,n
   real(8) :: elj
   real(8) :: r,sig,eps,qq

   n = size(pairs,1)
   
   elj = 0d0   
   do i = 1, 2
      do j = 3, n
         r = pairs(i,j)%rij
         sig = pairs(i,j)%ljSig
         eps = pairs(i,j)%ljEps
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
      end do
   end do
end function get_complexsolvent_lj_energy

subroutine get_complexsolvent_energy_with_H(pairs,elj,eel,e)
implicit none

   integer :: i,j,n
   real(8),intent(out) :: elj,eel,e
   real(8) :: r,sig,eps,qq
   type(AtomPairData),dimension(:,:),intent(in) :: pairs

   n = size(pairs,1)
   
   elj = 0d0   
   do i = 1, 2
      do j = 4, n
         r = pairs(i,j)%rij
         sig = pairs(i,j)%ljSig
         eps = pairs(i,j)%ljEps
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
      end do
   end do
   
   eel = 0d0   
   do i = 1, 3
      do j = 4, n
         r = pairs(i,j)%rij
         qq = pairs(i,j)%qq
         eel = eel + kCoulomb*qq/r
      end do
   end do
   
   e = elj + eel   
end subroutine get_complexsolvent_energy_with_H

subroutine get_solventsolvent_energy(pairs,elj,eel,eb,et) 
implicit none
   integer :: i,j,n
   real(8) :: r,sig,eps,qq
   real(8),intent(out) :: elj,eel,eb,et
   type(AtomPairData),dimension(:,:),intent(in) :: pairs

   n = size(pairs,1)
   
   eb = 0d0
   elj = 0d0
   eel = 0d0
   do i = 3, n-1, 2
      !eb = eb + 150d0*(pairs(i,i+1)%rij - 1.781d0)**2
      do j = i+2, n
         r = pairs(i,j)%rij
         sig = pairs(i,j)%ljSig
         eps = pairs(i,j)%ljEps
         qq = pairs(i,j)%qq
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
         eel = eel + kCoulomb*qq/r
         
         r = pairs(i+1,j)%rij
         sig = pairs(i+1,j)%ljSig
         eps = pairs(i+1,j)%ljEps
         qq = pairs(i+1,j)%qq
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
         eel = eel + kCoulomb*qq/r
      end do
   end do
   
   et = eb + elj + eel
end subroutine get_solventsolvent_energy

subroutine get_solventsolvent_energy_with_H(pairs,elj,eel,eb,et) 
implicit none
   integer :: i,j,n
   real(8) :: r,sig,eps,qq
   real(8),intent(out) :: elj,eel,eb,et
   type(AtomPairData),dimension(:,:),intent(in) :: pairs

   n = size(pairs,1)
   
   eb = 0d0
   elj = 0d0
   eel = 0d0
   do i = 4, n-1, 2
      !eb = eb + 150d0*(pairs(i,i+1)%rij - 1.781d0)**2
      do j = i+2, n
         r = pairs(i,j)%rij
         sig = pairs(i,j)%ljSig
         eps = pairs(i,j)%ljEps
         qq = pairs(i,j)%qq
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
         eel = eel + kCoulomb*qq/r
         
         r = pairs(i+1,j)%rij
         sig = pairs(i+1,j)%ljSig
         eps = pairs(i+1,j)%ljEps
         qq = pairs(i+1,j)%qq
         elj = elj + 4d0*eps*(sig**12/r**12-sig**6/r**6)
         eel = eel + kCoulomb*qq/r
      end do
   end do
   
   et = eb + elj + eel
end subroutine get_solventsolvent_energy_with_H

subroutine get_total_potential_energy(pairs,ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et)
implicit none
   
   real(8),intent(out) :: ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   
   ec = get_complex_energy(1d0,pairs(1,2)%rij)
   call get_complexsolvent_energy(pairs,ecslj,ecsel,ecst)
   call get_solventsolvent_energy(pairs,esslj,essel,essb,esst)
   et = ec + ecst + esst
end subroutine get_total_potential_energy

subroutine get_total_potential_energy_pbme(pairs,p,ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et)
use maproutines
implicit none
   
   real(8),intent(out) :: ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   type(QuantumStateData),intent(in) :: p
   
   integer :: i,nm
   real(8) :: trace
   real(8),dimension(:,:),allocatable :: vcsel,vcatt

   nm = size(p%eigenvalues)
   
   allocate(vcsel(1:nm,1:nm))
   allocate(vcatt(1:nm,1:nm))
   vcsel = 0d0
   vcatt = 0d0

   ec = get_complex_energy_repulsion_part(pairs(1,2)%rij)
   call make_matrix_traceless(p%hs,trace,vcatt)
   ec = ec + get_map_contribution(vcatt,p%mapFactor) + trace

   ecslj = get_complexsolvent_lj_energy(pairs)
   call make_matrix_traceless((p%vas+p%vbs+p%vhs),trace,vcsel)
   ecsel = get_map_contribution(vcsel,p%mapFactor) + trace
   ecst = ecslj + ecsel

   call get_solventsolvent_energy(pairs,esslj,essel,essb,esst)
   et = ec + ecst + esst
end subroutine get_total_potential_energy_pbme

subroutine get_total_potential_energy_fbts(pairs,p,ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et)
use maproutines
implicit none
   
   real(8),intent(out) :: ec,ecslj,ecsel,ecst,esslj,essel,essb,esst,et
   type(AtomPairData),dimension(:,:),intent(in) :: pairs
   type(QuantumStateData),intent(in) :: p
   
   integer :: i,nm
   real(8) :: trace
   real(8),dimension(:,:),allocatable :: vcsel,vcatt

   nm = size(p%eigenvalues)
   
   allocate(vcsel(1:nm,1:nm))
   allocate(vcatt(1:nm,1:nm))
   vcsel = 0d0
   vcatt = 0d0

   ec = get_complex_energy_repulsion_part(pairs(1,2)%rij)
   call make_matrix_traceless(p%hs,trace,vcatt)
   ec = ec + 0.5d0*(get_map_contribution(vcatt,p%mapFactor1) + &
                     get_map_contribution(vcatt,p%mapFactor2)) + trace

   ecslj = get_complexsolvent_lj_energy(pairs)
   call make_matrix_traceless((p%vas+p%vbs+p%vhs),trace,vcsel)
   ecsel = 0.5d0*(get_map_contribution(vcsel,p%mapFactor1) + &
                  get_map_contribution(vcsel,p%mapFactor2)) + trace
   ecst = ecslj + ecsel

   call get_solventsolvent_energy(pairs,esslj,essel,essb,esst)
   et = ec + ecst + esst
end subroutine get_total_potential_energy_fbts

end module energycalculation
