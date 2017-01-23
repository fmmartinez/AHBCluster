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
   
   m = nPointsGrid
   pp = pp + f%gridPointValue(m)*g%gridPointValue(m)*h%gridPointValue(m)

   do i = 2, nPointsGrid-1
      pp = pp + f%gridPointValue(i)*g%gridPointValue(i)*h%gridPointValue(i)
   end do

   pp = pp*binWidth

end function

subroutine get_overlap_matrix(phi,S)
implicit none
   real(8),dimension(:,:),intent(out) :: S
   type(BasisFunction),dimension(:),intent(in) :: phi
   
   integer :: i,j,n
   type(EvalOnGridFunction) :: valueOfOne

   valueofOne%gridPointValue(1:nPointsGrid) = 1d0

   n = size(phi)

   do i = 1, n
      do j = 1, n
         S(i,j) = integrate_trapezoid_rule(phi(i),valueOfOne,phi(j))
      end do
   end do
end subroutine get_overlap_matrix

subroutine get_phi_KineticEnergy_phi_matrix(phi,d2p,K)
implicit none
   real(8),dimension(:,:),intent(out) :: K 
   type(BasisFunction),dimension(:),intent(in) :: phi,d2p
   
   integer :: i,j,n
   type(EvalOnGridFunction) :: valueOfOne

   valueofOne%gridPointValue(1:nPointsGrid) = 1d0

   n = size(phi)
   
   do i = 1, n
      do j = 1, n
         K(i,j) = -0.0479d0*integrate_trapezoid_rule(d2p(i),valueOfOne,phi(j))
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

   do i = 1, nPointsGrid
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

subroutine get_double_derivative_basis_functions_on_each_well(cov,ion)
implicit none
   type(BasisFunction),dimension(:),intent(inout) :: cov,ion
   
   integer :: i,j,nbc,nbi
   real(8) :: q
   
   nbc = size(cov)
   nbi = size(ion)

   do i = 1, nPointsGrid
      q = lowerLimit + (i-1)*binWidth
      !j enters calculation as j-1, because the first index (j=1) corresponds to 0, 
      !which by definition is ground state
      do j = 1, nbc
         cov(j)%gridPointValue(i) = eval_d2_harmonic_oscillator_wavefunction(j-1,q-covMinWell)
      end do
      do j = 1, nbi
         ion(j)%gridPointValue(i) = eval_d2_harmonic_oscillator_wavefunction(j-1,q-ionMinWell)
      end do
   end do
end subroutine get_double_derivative_basis_functions_on_each_well

subroutine initialize_basis_functions_on_each_well(cov,ion)
implicit none
   type(BasisFunction),dimension(:),intent(inout) :: cov,ion
   
   integer :: i,j,nbc,nbi
   real(8) :: q
   
   nbc = size(cov)
   nbi = size(ion)

   do i = 1, nPointsGrid
      q = lowerLimit + (i-1)*binWidth
      !j enters calculation as j-1, because the first index (j=1) corresponds to 0, 
      !which by definition is ground state
      do j = 1, nbc
         cov(j)%gridPointValue(i) = eval_harmonic_oscillator_wavefunction(j-1,q-covMinWell)
      end do
      do j = 1, nbi
         ion(j)%gridPointValue(i) = eval_harmonic_oscillator_wavefunction(j-1,q-ionMinWell)
      end do
   end do
end subroutine initialize_basis_functions_on_each_well

function eval_harmonic_oscillator_wavefunction(i,dx) result(phi)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx

   real(8) :: phi
   real(8) :: hermite,factorial,exponent1

   hermite = eval_hermite_polynomial(i,alpha*dx)
   factorial = eval_factorial(i)
   exponent1 = exp(-1d0*(alpha*dx)**2/2d0)

   phi = sqrt(alpha/(2d0**i*factorial*pisqrt))*hermite*exponent1

end function eval_harmonic_oscillator_wavefunction

function eval_d2_harmonic_oscillator_wavefunction(i,dx) result(d2p)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx

   real(8) :: d2p,a,f,g,df,dg,d2f,d2g,factorial

   call eval_derivatives_gaussian_function(dx,g,dg,d2g)
   call eval_derivativez_hermite_polynomial(i,alpha*dx,f,df,d2f)
   factorial = eval_factorial(i)

   a = sqrt(alpha/(2d0**i*factorial*pisqrt))
   d2p = a*f*d2g + 2d0*a*dg*df + a*g*d2f

end function eval_d2_harmonic_oscillator_wavefunction

subroutine eval_derivatives_gaussian_function(dx,d0,d1,d2)
!the form of the gaussian is exp(-dx^2*c^2/2)
implicit none
   real(8),intent(in) :: dx
   real(8),intent(out) :: d0,d1,d2
   
   d0 = exp(-0.5d0*alpha**2*dx**2)
   
   d1 = -d0*alpha**2*dx

   d2 = d0*alpha**2*(alpha**2*dx**2 - 1d0)

end subroutine eval_derivatives_gaussian_function

subroutine eval_derivativez_hermite_polynomial(i,dx,d0,d1,d2)
implicit none
   integer,intent(in) :: i
   real(8),intent(in) :: dx
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
   
   e(1:m) = w(1:m)
   l(1:n,1:m) = Htemp(1:n,1:m)

   deallocate(Stemp)
   deallocate(Htemp)
end subroutine

end module quantumcalculations
