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
end subroutine

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
      print *, 'warning. Hermite polynomial with negative index'
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

end module quantumcalculations
