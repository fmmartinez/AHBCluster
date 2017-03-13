module maproutines
use definitions

contains

subroutine do_mapping_variables_sampling_box_muller(stream,p)
use mkl_vsl_type
use mkl_vsl
implicit none
   type(vsl_stream_state),intent(in) :: stream
   type(QuantumStateData),intent(inout) :: p
   
   integer :: nMap,errcode
   
   nMap = size(p%rm)

   errcode = vdrnggaussian(method,stream,nMap,p%rm,0d0,hbar)!sqrt(hbar/2d0))
   errcode = vdrnggaussian(method,stream,nMap,p%pm,0d0,hbar)!sqrt(hbar/2d0))

end subroutine do_mapping_variables_sampling_box_muller

subroutine do_mapping_variables_sampling(stream,p,state0)
use mkl_vsl_type
use mkl_vsl
implicit none
   integer,intent(in) :: state0
   type(vsl_stream_state),intent(in) :: stream
   type(QuantumStateData),intent(inout) :: p
   
   integer :: nMap,errcode,i
   integer,parameter :: method1 = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
   real(8) :: limit
   real(8),dimension(:),allocatable :: angles

   nMap = size(p%rm)

   allocate(angles(1:nMap))

   errcode = vdrnguniform(method1,stream,nMap,angles,0d0,2d0*pi)

   if (nMap == 1) then
      p%rm(1) = sqrt(3d0*hbar)*cos(angles(1))
      p%pm(1) = sqrt(3d0*hbar)*sin(angles(1))
   else
      do i = 1, nMap
         if (i == state0) then
            p%rm(i) = sqrt(3d0*hbar)*cos(angles(i))
            p%pm(i) = sqrt(3d0*hbar)*sin(angles(i))
         else
            p%rm(i) = sqrt(hbar)*cos(angles(i))
            p%pm(i) = sqrt(hbar)*sin(angles(i))
         end if
      end do
   end if

   
   !limit = sqrt(hbar)

   !errcode = vdrnguniform(method1,stream,nMap,p%rm,-limit,limit)

   !if (nMap == 1) then
   !   p%rm(1) = sqrt(3d0)*p%rm(1)
   !   p%pm(1) = sqrt(3d0*hbar - p%rm(1)**2)
   !else
   !   do i = 1, nMap
   !      if (i == state0) then
   !         p%rm(i) = sqrt(3d0)*p%rm(i)
   !         p%pm(i) = sqrt(3d0*hbar - p%rm(i)**2)
   !      else
   !         !rm stays same
   !         p%pm(i) = sqrt(hbar - p%rm(i)**2)
   !      end if
   !   end do
   !end if

end subroutine do_mapping_variables_sampling

subroutine do_coherent_state_variables_sampling(stream,p)
use mkl_vsl_type
use mkl_vsl
implicit none
   type(vsl_stream_state),intent(in) :: stream
   type(QuantumStateData),intent(inout) :: p
   
   integer :: nMap,errcode
   
   nMap = size(p%rm)

   errcode = vdrnggaussian(method,stream,nMap,p%p1,0d0,hbar)
   errcode = vdrnggaussian(method,stream,nMap,p%q2,0d0,hbar)
   
   errcode = vdrnggaussian(method,stream,nMap,p%p2,0d0,hbar)
   errcode = vdrnggaussian(method,stream,nMap,p%q2,0d0,hbar)

end subroutine do_coherent_state_variables_sampling

function get_map_contribution(mat,mapFactor) result(y)
implicit none
   real(8),dimension(:,:),intent(in) :: mapFactor,mat
   
   integer :: i,j,nm
   real(8) :: y
   
   nm = size(mapFactor,1)

   y = 0d0
   do i = 1, nm
      y = y + mat(i,i)*mapFactor(i,i)
   end do

   do i = 1, nm-1
      do j = i+1, nm
         y = y + mat(i,j)*mapFactor(i,j)*2d0
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
      pbme%mapFactor(i,i) = 0.5d0*(pbme%rm(i)**2 + pbme%pm(i)**2 - hbar)/hbar
   end do
   
   do i = 1, nm-1
      do j = i+1, nm
         pbme%mapFactor(i,j) = 0.5d0*(pbme%rm(i)*pbme%rm(j) + pbme%pm(i)*pbme%pm(j))/hbar
         pbme%mapFactor(j,i) = pbme%mapFactor(i,j)
      end do
   end do
end subroutine get_mapFactor

subroutine get_mapFactor_traceless(pbme)
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

end subroutine get_mapFactor_traceless

subroutine get_covarFactor(fbts)
implicit none
   type(QuantumStateData),intent(inout) :: fbts
   integer :: i,j,nm

   nm = size(fbts%q1)

   fbts%mapFactor1 = 0d0
   fbts%mapFactor2 = 0d0
   do i = 1, nm
      do j = 1, nm
         fbts%mapFactor1(i,j) = 0.5d0*(fbts%q1(i)*fbts%q1(j) + fbts%p1(i)*fbts%p1(j))/hbar
         fbts%mapFactor2(i,j) = 0.5d0*(fbts%q2(i)*fbts%q2(j) + fbts%p2(i)*fbts%p2(j))/hbar
      end do
   end do

end subroutine get_covarFactor

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

function get_apparent_rAH(pbme) result(aqA)
use quantumcalculations
implicit none
   type(QuantumStateData),intent(in) :: pbme
   
   integer :: i,j,nb,nm
   real(8) :: aqA, trace
   real(8),dimension(:),allocatable :: li,lj
   real(8),dimension(:,:),allocatable :: q,qt
   
   nb = size(pbme%lambda,1)
   nm = size(pbme%rm)
   
   allocate(li(1:nb))
   allocate(lj(1:nb))
   allocate(q(1:nm,1:nm))
   allocate(qt(1:nm,1:nm))
   
   do i = 1, nm
      do j = 1, nm
         li = pbme%lambda(1:nb,i)
         lj = pbme%lambda(1:nb,j)
         q(i,j) = get_lambda_f_lambda_matrix_element(li,lj,pbme%prAHp)
      end do
   end do

   call make_matrix_traceless(q,trace,qt)
   aqA = get_map_contribution(qt,pbme%mapFactor) + trace
end function get_apparent_rAH

end module maproutines
