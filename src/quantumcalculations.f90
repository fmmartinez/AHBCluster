module quantumcalculations
use definitions
implicit none

contains

function integrate_trapezoid_rule(f,g,h) result(pp)
!functions evaluates \int f_i g h_i
implicit none
   type(BasisFunction),intent(in) :: f, h
   type(EvalOnGridFunction),intent(in) :: g
   
   integer :: i,m
   real(8) :: pp

   pp = f%gridPointValue(1)*g%gridPointValue(1)*h%gridPointValue(1)
   
   m = nPointsGrid+1
   pp = pp + f%gridPointValue(m)*g%gridPointValue(m)*h%gridPointValue(m)
   
   pp = 0.5d0*pp

   do i = 2, nPointsGrid
      pp = pp + f%gridPointValue(i)*g%gridPointValue(i)*h%gridPointValue(i)
   end do

   pp = pp*binWidth

end function

!---matrix elements and matrix assembly
!--- f stands for general function not force
!--- <phi | f | phi> matrix elements
subroutine get_overlap_matrix(phi,S)
implicit none
   real(8),dimension(:,:),intent(out) :: S
   type(BasisFunction),dimension(:),intent(in) :: phi
   
   integer :: i,j,n
   type(EvalOnGridFunction) :: valueOfOne

   valueofOne%gridPointValue(1:nPointsGrid+1) = 1d0

   n = size(phi)

   do i = 1, n
      do j = 1, n
         S(i,j) = integrate_trapezoid_rule(phi(i),valueOfOne,phi(j))
      end do
   end do
end subroutine get_overlap_matrix

!---These are used to construct the subsystem and diagonalize
subroutine get_phi_KineticEnergy_phi_matrix(phi,d2p,K)
implicit none
   real(8),dimension(:,:),intent(out) :: K 
   type(BasisFunction),dimension(:),intent(in) :: phi,d2p
   
   integer :: i,j,n
   type(EvalOnGridFunction) :: valueOfOne

   valueofOne%gridPointValue(1:nPointsGrid+1) = 1d0

   n = size(phi)
   
   do i = 1, n
      do j = 1, n
         K(i,j) = -0.0479d0*integrate_trapezoid_rule(phi(i),valueOfOne,d2p(j))
      end do
   end do
end subroutine get_phi_KineticEnergy_phi_matrix

subroutine get_phi_Vsubsystem_phi_matrix(phi,rab,V)
use energycalculation
implicit none
   real(8),intent(in) :: rab
   real(8),dimension(:,:),intent(out) :: V
   type(basisFunction),dimension(:),intent(in) :: phi
   
   integer :: i,j,n
   real(8) :: q
   type(EvalOnGridFunction) :: vh

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      vh%gridPointValue(i) = get_complex_energy_attraction_part(q,rab)
   end do

   n = size(phi)

   do i = 1, n
      do j = 1, n
         V(i,j) = integrate_trapezoid_rule(phi(i),vh,phi(j))
      end do
   end do

end subroutine get_phi_Vsubsystem_phi_matrix

subroutine get_phi_rAH_phi_matrix(p)
implicit none
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,n
   type(EvalOnGridFunction) :: rHS
   
   n = size(p%phi)

   do i = 1, nPointsGrid+1
      rHS%gridPointValue(i) = p%gridHSolvent(1)%gridPoint(i)%rij
   end do
   
   do i = 1, n
      do j = 1, n
         p%prAHp(i,j) = integrate_trapezoid_rule(p%phi(i),rHS,p%phi(j))
      end do
   end do
end subroutine get_phi_rAH_phi_matrix

!---various matrix elements used through energy and force calculations
subroutine get_phi_charge_AB_phi_matrix(p)
implicit none
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,n
   real(8) :: q,fr
   type(EvalOnGridFunction) :: qA,qB

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      fr = 0.5d0*(1d0 + (q - 1.43d0)/sqrt((q-1.43d0)**2+0.125d0**2))
      qA%gridPointValue(i) = (1d0-fr)*(-0.5d0)+fr*(-1d0)
      qB%gridPointValue(i) = fr*(0.5d0)
   end do

   n = size(p%phi)

   do i = 1, n
      do j = 1, n
         p%pqAp(i,j) = integrate_trapezoid_rule(p%phi(i),qA,p%phi(j))
         p%pqBp(i,j) = integrate_trapezoid_rule(p%phi(i),qB,p%phi(j))
      end do
   end do
end subroutine get_phi_charge_AB_phi_matrix

subroutine get_phi_d_VAH_phi_matrix(p)
implicit none
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,n
   real(8) :: q,rah
   type(EvalOnGridFunction) :: vh

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      rah = q
      vh%gridPointValue(i) = &
         -d*exp(-na*(rah-da)**2/(2d0*rah))*(na*(rah-da)/rah)*((rah-da)/(2d0*rah)-1d0)
   end do

   n = size(p%phi)

   do i = 1, n
      do j = 1, n
         p%pAHp(i,j) = integrate_trapezoid_rule(p%phi(i),vh,p%phi(j))
      end do
   end do

end subroutine get_phi_d_VAH_phi_matrix

subroutine get_phi_d_VBH_phi_matrix(p,rab)
implicit none
   real(8),intent(in) :: rab
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,n
   real(8) :: q,rbh
   type(EvalOnGridFunction) :: vh

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      rbh = rab - q
      vh%gridPointValue(i) = &
         -c*d*exp(-nb*(rbh-db)**2/(2d0*rbh))*(nb*(rbh-db)/rbh)*((rbh-db)/(2d0*rbh)-1d0)
   end do
!print *, '--',rbh,vh%gridPointValue(nPointsGrid+1)
   n = size(p%phi)
   
   do i = 1, n
      do j = 1, n
         p%pBHp(i,j) = integrate_trapezoid_rule(p%phi(i),vh,p%phi(j))
      end do
   end do

end subroutine get_phi_d_VBH_phi_matrix

subroutine get_phi_inv_r_HS_phi_matrix(p)
implicit none
   type(QuantumStateData),intent(inout) :: p

   integer :: i,j,k,n,nAtoms
   type(EvalOnGridFunction) :: inverse_rHS
   
   n = size(p%phi)
   nAtoms = size(p%gridHSolvent)

   do k = 1, nAtoms
      do i = 1, nPointsGrid+1
         inverse_rHS%gridPointValue(i) = 1d0/p%gridHSolvent(k)%gridPoint(i)%rij
      end do
      
      do i = 1, n
         do j = 1, n
            p%pirp(k)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse_rHS,p%phi(j))
         end do
      end do
   end do
end subroutine get_phi_inv_r_HS_phi_matrix

!subroutine get_phi_inv_r2_HS_phi_matrix(p)
!!inv2 is a vector
!implicit none
!   type(QuantumStateData),intent(inout) :: p
!   
!   integer :: i,j,k,n,nAtoms
!   type(EvalOnGridFunction),dimension(1:3) :: inverse_rHS_2
!   
!   n = size(p%phi)
!   nAtoms = size(p%gridHSolvent)
!
!   do k = 1, nAtoms
!      do i = 1, nPointsGrid+1
!         !negative make pointing from grid to atom, negative from derivative
!         !makes it positive
!         inverse_rHS_2%gridPointValue(i) = &
!            p%gridHSolvent(k)%gridPoint(i)%vectorij/(p%gridHSolvent(k)%gridPoint(i)%rij**3)
!      end do
!      
!      do i = 1, n
!         do j = 1, n
!            p%pir2p(k,1)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse_rHS_2(1),p%phi(j))
!            p%pir2p(k,2)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse_rHS_2(2),p%phi(j))
!            p%pir2p(k,3)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse_rHS_2(3),p%phi(j))
!         end do
!      end do
!   end do
!end subroutine get_phi_inv_r2_HS_phi_matrix

subroutine get_phi_inv_r3_HS_phi_matrix(p)
!inv2 is a vector
implicit none
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,n,nAtoms
   type(EvalOnGridFunction) :: inverse3
   
   n = size(p%phi)
   nAtoms = size(p%gridHSolvent)

   do k = 1, nAtoms
      do i = 1, nPointsGrid+1
         inverse3%gridPointValue(i) = 1d0/(p%gridHSolvent(k)%gridPoint(i)%rij**3)
      end do

      do i = 1, n
         do j = 1, n
            p%pir3p(k)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse3,p%phi(j))
         end do
      end do
   end do
end subroutine get_phi_inv_r3_HS_phi_matrix

subroutine get_phi_rc_inv_r3_HS_phi_matrix(rab,p)
!inv2 is a vector
implicit none
   real(8),intent(in) :: rab
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,n,nAtoms
   real(8) :: rch,rsh
   type(EvalOnGridFunction) :: inverse3
   
   n = size(p%phi)
   nAtoms = size(p%gridHSolvent)

   do k = 1, nAtoms
      do i = 1, nPointsGrid+1
         rch = rab*(0.3882d0) - p%gridHSolvent(1)%gridPoint(i)%rij
         rsh = p%gridHSolvent(k)%gridPoint(i)%rij
         inverse3%gridPointValue(i) = rch/(rsh**3)
      end do
      
      do i = 1, n
         do j = 1, n
            p%pcr3p(k)%mat(i,j) = integrate_trapezoid_rule(p%phi(i),inverse3,p%phi(j))
         end do
      end do
   end do
end subroutine get_phi_rc_inv_r3_HS_phi_matrix

!--- <lambda | f | lambda > matrix elements
subroutine get_lambda_h_lambda_matrix(at,pair,p)
implicit none
   type(Atom),dimension(:),intent(in) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(QuantumStateData),intent(inout) :: p

   integer :: i,m

   m = size(p%eigenvalues)
   
   p%h = 0d0   
   call get_lambda_hs_lambda_matrix(p)
   !do i = 1, m
   !   p%h(i,i) = p%eigenvalues(i)
   !end do
   p%h = p%hs
   
   call get_lambda_VASol_lambda_matrix(at,pair,p)
   call get_lambda_VBSol_lambda_matrix(at,pair,p)
   call get_lambda_VHSol_lambda_matrix(at,p)
   
   p%h = p%h + p%vas + p%vbs + p%vhs

end subroutine get_lambda_h_lambda_matrix

subroutine get_lambda_hs_lambda_matrix(p)
implicit none
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,nb,nm
   real(8) :: vc
   real(8),dimension(:),allocatable :: li,lj
   real(8),dimension(:,:),allocatable :: phihsphi
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))
   allocate(phihsphi(1:nb,1:nb))

   p%hs = 0d0
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         phihsphi = p%phiKphi + p%phiVsphi
         vc = get_lambda_f_lambda_matrix_element(li,lj,phihsphi)
         p%hs(i,j) = p%hs(i,j) + vc
      end do
   end do
   
   deallocate(phihsphi)
   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_hs_lambda_matrix

subroutine get_lambda_VHSol_lambda_matrix(at,p)
implicit none
   real(8),parameter :: hCharge = 0.5d0
   type(Atom),dimension(:),intent(in) :: at
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,na,nb,nm
   real(8) :: vc,prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   na = size(at)
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   p%vhs = 0d0
   do k = 3, na
      prefactor = kCoulomb*at(k)%charge*hCharge
      do i = 1, nm
         do j = 1, nm
            li = p%lambda(1:nb,i)
            lj = p%lambda(1:nb,j)
            vc = get_lambda_f_lambda_matrix_element(li,lj,p%pirp(k)%mat)
            p%vhs(i,j) = p%vhs(i,j) + prefactor*vc
         end do
      end do
   end do
   
   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_VHSol_lambda_matrix

subroutine get_lambda_VASol_lambda_matrix(at,pair,p)
implicit none
   type(Atom),dimension(:),intent(in) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,na,nb,nm
   real(8) :: vc,prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   na = size(at)
   nb = size(p%lambda,1)
   nm = size(p%rm)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   p%vas = 0d0
   do k = 3, na
      prefactor = (kCoulomb*at(k)%charge/pair(1,k)%rij)
      do i = 1, nm
         do j = 1, nm
            li = p%lambda(1:nb,i)
            lj = p%lambda(1:nb,j)
            vc = get_lambda_f_lambda_matrix_element(li,lj,p%pqAp)
            p%vas(i,j) = p%vas(i,j) + prefactor*vc
         end do
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_VASol_lambda_matrix

subroutine get_lambda_VBSol_lambda_matrix(at,pair,p)
implicit none
   type(Atom),dimension(:),intent(in) :: at
   type(AtomPairData),dimension(:,:),intent(in) :: pair
   type(QuantumStateData),intent(inout) :: p
   
   integer :: i,j,k,na,nb,nm
   real(8) :: vc,prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   na = size(at)
   nb = size(p%lambda,1)
   nm = size(p%rm)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   p%vbs = 0d0
   do k = 3, na
      prefactor = (kCoulomb*at(k)%charge/pair(2,k)%rij)
      do i = 1, nm
         do j = 1, nm
            li = p%lambda(1:nb,i)
            lj = p%lambda(1:nb,j)
            vc = get_lambda_f_lambda_matrix_element(li,lj,p%pqBp)
            p%vbs(i,j) = p%vbs(i,j) + prefactor*vc
         end do
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_VBSol_lambda_matrix

subroutine get_lambda_q_lambda_matrix(lambda,phifphi,v)
implicit none
   real(8),dimension(:,:),intent(in) :: lambda, phifphi
   real(8),dimension(:,:),intent(out) :: v
   
   integer :: i,j,b,nb,nm
   real(8) :: vc
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(lambda,1)
   nm = size(v,1)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   v = 0d0
   do i = 1, nm
      do j = 1, nm
         li = lambda(1:nb,i)
         lj = lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,phifphi)
         v(i,j) = v(i,j) + vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_q_lambda_matrix

!--- df stands for derivative with respect to R of f
!---<lambda | df | lambda > matrix elements
subroutine get_lambda_d_VASol_lambda_matrix(at,pair,p,dv)
implicit none
   type(Atom),intent(in) :: at
   type(AtomPairData),intent(in) :: pair
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc, prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   prefactor = -kCoulomb*at%charge/pair%rij**2
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pqAp)
         dv(i,j) = dv(i,j) + prefactor*vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VASol_lambda_matrix

subroutine get_lambda_d_VBSol_lambda_matrix(at,pair,p,dv)
implicit none
   type(Atom),intent(in) :: at
   type(AtomPairData),intent(in) :: pair
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc, prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   prefactor = -kCoulomb*at%charge/pair%rij**2
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pqBp)
         dv(i,j) = dv(i,j) + prefactor*vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VBSol_lambda_matrix

subroutine get_lambda_d_VHSol_lambda_matrix(at,lambda,phifphi,dvx,dvy,dvz)
!phifphi is a generic name, in this case stored in that variable is
!<phi| rvec/r^3 |phi>
implicit none
   real(8),parameter :: hCharge = 0.5d0
   type(Atom),intent(in) :: at
   type(MatrixList),dimension(:),intent(in) :: phifphi
   real(8),dimension(:,:),intent(in) :: lambda
   real(8),dimension(:,:),intent(out) :: dvx,dvy,dvz
   
   integer :: i,j,nb,nm
   real(8) :: vc,prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(lambda,1)
   nm = size(dvx,1)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dvx = 0d0
   dvy = 0d0
   dvz = 0d0
   prefactor = kCoulomb*at%charge*hCharge
   do i = 1, nm
      do j = 1, nm
         li = lambda(1:nb,i)
         lj = lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,phifphi(1)%mat)
         dvx(i,j) = dvx(i,j) + prefactor*vc
         vc = get_lambda_f_lambda_matrix_element(li,lj,phifphi(2)%mat)
         dvy(i,j) = dvy(i,j) + prefactor*vc
         vc = get_lambda_f_lambda_matrix_element(li,lj,phifphi(3)%mat)
         dvz(i,j) = dvz(i,j) + prefactor*vc
      end do
   end do
   
   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VHSol_lambda_matrix

subroutine get_lambda_d_VCoMSol_lambda_matrix(k,at,p,dv)
implicit none
   real(8),parameter :: chargeH = 0.5d0
   
   integer,intent(in) :: k
   type(Atom),intent(in) :: at
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc, prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   prefactor = -kCoulomb*at%charge*chargeH
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pir3p(k)%mat)
         dv(i,j) = dv(i,j) + prefactor*vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VCoMSol_lambda_matrix

subroutine get_lambda_d_VCoMH_lambda_matrix(k,at,p,dv)
implicit none
   real(8),parameter :: chargeH = 0.5d0
   
   integer,intent(in) :: k
   type(Atom),intent(in) :: at
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc, prefactor
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   prefactor = -kCoulomb*at%charge*chargeH
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pcr3p(k)%mat)
         dv(i,j) = dv(i,j) + prefactor*vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VCoMH_lambda_matrix

subroutine get_lambda_d_VAH_lambda_matrix(p,dv)
implicit none
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc 
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pAHp)
         dv(i,j) = dv(i,j) + vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VAH_lambda_matrix

subroutine get_lambda_d_VBH_lambda_matrix(p,dv)
implicit none
   type(QuantumStateData),intent(in) :: p
   real(8),dimension(:,:),intent(out) :: dv
   
   integer :: i,j,nb,nm
   real(8) :: vc 
   real(8),dimension(:),allocatable :: li,lj
   
   nb = size(p%lambda,1)
   nm = size(p%eigenvalues)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))

   dv = 0d0
   do i = 1, nm
      do j = 1, nm
         li = p%lambda(1:nb,i)
         lj = p%lambda(1:nb,j)
         vc = get_lambda_f_lambda_matrix_element(li,lj,p%pBHp)
         dv(i,j) = dv(i,j) + vc
      end do
   end do

   deallocate(lj)
   deallocate(li)
end subroutine get_lambda_d_VBH_lambda_matrix

!---individual  <lambda | f | lambda > matrix elements calculator
function get_lambda_f_lambda_matrix_element(l1,l2,phifphi) result(f)
implicit none
   real(8),dimension(:),intent(in) :: l1,l2
   real(8),dimension(:,:),intent(in) :: phifphi
   
   integer :: nb,i,j
   real(8) :: f

   nb = size(l1)
   
   f = 0d0
   do i = 1, nb
      do j = 1, nb
         f = f + l1(i)*phifphi(i,j)*l2(j)
      end do
   end do
end function get_lambda_f_lambda_matrix_element

!---basis functions preparation
subroutine get_double_derivative_basis_functions_on_each_well(cov,ion)
implicit none
   type(BasisFunction),dimension(:),intent(inout) :: cov,ion
   
   integer :: i,j,nbc,nbi
   real(8) :: q,alpha
   
   nbc = size(cov)
   nbi = size(ion)

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      !j enters calculation as j-1, because the first index (j=1) corresponds to 0, 
      !which by definition is ground state
      do j = 1, nbc
         alpha = alphaCov
         cov(j)%gridPointValue(i) = eval_d2_harmonic_oscillator_wavefunction(j-1,q-covMinWell,alpha)
      end do
      do j = 1, nbi
         alpha = alphaIon
         ion(j)%gridPointValue(i) = eval_d2_harmonic_oscillator_wavefunction(j-1,q-ionMinWell,alpha)
      end do
   end do
end subroutine get_double_derivative_basis_functions_on_each_well

subroutine initialize_basis_functions_on_each_well(cov,ion)
implicit none
   type(BasisFunction),dimension(:),intent(inout) :: cov,ion
   
   integer :: i,j,nbc,nbi
   real(8) :: q,alpha
   
   nbc = size(cov)
   nbi = size(ion)

   do i = 1, nPointsGrid+1
      q = lowerLimit + (i-1)*binWidth
      !j enters calculation as j-1, because the first index (j=1) corresponds to 0, 
      !which by definition is ground state
      do j = 1, nbc
         alpha = alphaCov
         cov(j)%gridPointValue(i) = eval_harmonic_oscillator_wavefunction(j-1,q-covMinWell,alpha)
      end do
      do j = 1, nbi
         alpha = alphaIon
         ion(j)%gridPointValue(i) = eval_harmonic_oscillator_wavefunction(j-1,q-ionMinWell,alpha)
      end do
   end do
end subroutine initialize_basis_functions_on_each_well

!---elementary calculations
function eval_harmonic_oscillator_wavefunction(i,dx,alpha) result(phi)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx,alpha

   real(8) :: phi
   real(8) :: hermite,factorial,exponent1

   hermite = eval_hermite_polynomial(i,alpha*dx)
   factorial = eval_factorial(i)
   exponent1 = exp(-1d0*(alpha*dx)**2/2d0)

   phi = sqrt(alpha/(2d0**i*factorial*pisqrt))*hermite*exponent1

end function eval_harmonic_oscillator_wavefunction

function eval_d2_harmonic_oscillator_wavefunction(i,dx,alpha) result(d2p)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx,alpha

   real(8) :: d2p,a,f,g,df,dg,d2f,d2g,factorial

   call eval_derivatives_gaussian_function(alpha,dx,g,dg,d2g)
   call eval_derivativez_hermite_polynomial(i,alpha,alpha*dx,f,df,d2f)
   factorial = eval_factorial(i)

   a = sqrt(alpha/(2d0**i*factorial*pisqrt))
   d2p = a*f*d2g + 2d0*a*dg*df + a*g*d2f

end function eval_d2_harmonic_oscillator_wavefunction

subroutine eval_derivatives_gaussian_function(alpha,dx,d0,d1,d2)
!the form of the gaussian is exp(-dx^2*c^2/2)
implicit none
   real(8),intent(in) :: dx,alpha
   real(8),intent(out) :: d0,d1,d2
   
   d0 = exp(-0.5d0*alpha**2*dx**2)
   
   d1 = -d0*alpha**2*dx

   d2 = d0*alpha**2*(alpha**2*dx**2 - 1d0)

end subroutine eval_derivatives_gaussian_function

subroutine eval_derivativez_hermite_polynomial(i,alpha,dx,d0,d1,d2)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx,alpha
   real(8),intent(out) :: d0,d1,d2

   d0 = eval_hermite_polynomial(i,dx)

   d1 = 2d0*i*eval_hermite_polynomial(i-1,dx)*alpha

   d2 = 4d0*i*(i-1)*eval_hermite_polynomial(i-2,dx)*alpha**2

end subroutine eval_derivativez_hermite_polynomial

function eval_hermite_polynomial(degree,x) result(h)
implicit none
   integer,intent(in) :: degree
   real(8),intent(in) :: x
   
   integer :: i
   real(8) :: h,a,c,yn,y1,y0

   a = 2d0

   y0 = 1d0
   y1 = 2d0*x

   if (degree == 0) then
      h = y0
   else if (degree == 1) then
      h = y1
   else if (degree >= 2) then
      do i = 2, degree
         c = 2d0*(i - 1d0)
         yn = a*x*y1 - c*y0

         y0 = y1
         y1 = yn
      end do

      h = yn
   else if (degree < 0) then
      h = 0d0
   end if
end function eval_hermite_polynomial

function eval_factorial(n) result(f)
implicit none
   integer,intent(in) :: n

   real(8) :: f
   integer :: i

   f = 1
   do i = 1, n
      f = f*i
   end do
end function eval_factorial

!-- diagonalizations
subroutine get_subsystem_lambdas(H,S,l,e)
implicit none
   real(8),dimension(:),intent(out) :: e
   real(8),dimension(:,:),intent(in) :: H,S
   real(8),dimension(:,:),intent(out) :: l
   
   integer :: n,m,lwork,info_1
   real(8),dimension(:),allocatable :: w,work
   real(8),dimension(:,:),allocatable :: Htemp,Stemp
   
   n = size(H,1)
   m = size(e)

   allocate(Htemp(1:n,1:n))
   allocate(Stemp(1:n,1:n))
   allocate(w(1:n))
   allocate(work(1:n*50))
   
   lwork = n*40
   
   Htemp = H
   Stemp = S
   call dsygv(1,'V','U',n,Htemp,n,Stemp,n,w,work,lwork,info_1)
   
   e = w
   
   l = Htemp
   
   deallocate(Stemp)
   deallocate(Htemp)
end subroutine

end module quantumcalculations
