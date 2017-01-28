module maproutines
use definitions

contains

subroutine do_mapping_variables_sampling(stream,p)
use mkl_vsl_type
use mkl_vsl
implicit none
   type(vsl_stream_state),intent(inout) :: stream
   type(QuantumStateData),intent(inout) :: p
   
   integer :: nMap,errcode
   
   nMap = size(p%rm)

   errcode = vdrnggaussian(method,stream,nMap,p%rm,0d0,hbar)
   errcode = vdrnggaussian(method,stream,nMap,p%pm,0d0,hbar)

end subroutine do_mapping_variables_sampling

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

subroutine get_mapFactor(pbme)
implicit none
   type(QuantumStateData),intent(inout) :: pbme
   integer :: i,j,nm

   nm = size(pbme%rm)

   pbme%mapFactor = 0d0
   do i = 1, nm
      do j = 1, nm
         pbme%mapFactor(i,j) = 0.5d0*(pbme%rm(i)*pbme%rm(j) + pbme%pm(i)*pbme%pm(j))/hbar
      end do
   end do

end subroutine get_mapFactor

subroutine make_matrix_traceless(matrix,traceN,tracelessmatrix)
implicit none
   real(8),intent(out) :: traceN
   real(8),dimension(:,:),intent(in) :: matrix
   real(8),dimension(:,:),intent(out) :: tracelessMatrix
   
   integer :: i,n

   n = size(matrix,1)

   traceN = 0d0
   do i = 1, n
      traceN = traceN + matrix(i,i)
   end do

   tracelessMatrix = matrix
   do i = 1, n
      tracelessMatrix(i,i) = tracelessMatrix(i,i) - traceN/real(n)
   end do

   traceN = traceN/real(n)
end subroutine make_matrix_traceless

end module maproutines
