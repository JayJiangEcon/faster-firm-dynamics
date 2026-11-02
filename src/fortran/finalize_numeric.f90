! Local finishing preserves the Mac log library and Fortran reduction order.
! The model policy, distribution and simulation paths are supplied by CUDA.
module numerical_finishing
use iso_c_binding, only: c_double
implicit none
contains
subroutine direct_objective(moments,targets,weights,objective) bind(C,name='direct_objective')
 real(c_double), intent(in) :: moments(17),targets(17),weights(17,17)
 real(c_double), intent(out) :: objective
 real(c_double) :: momdiffvec(16),weightmat(16,16),quadformvec(16)
 integer :: i,j
 do i=1,16
  momdiffvec(i)=targets(i+1)-moments(i+1)
  do j=1,16
   weightmat(i,j)=weights(i+1,j+1)
  end do
 end do
 quadformvec=matmul(weightmat,momdiffvec)
 objective=dot_product(momdiffvec,quadformvec)
end subroutine
subroutine macro_variance(output,variance) bind(C,name='macro_variance')
 real(c_double), intent(in) :: output(302)
 real(c_double), intent(out) :: variance
 real(c_double) :: growth,mean,second
 integer :: t
 mean=0.0_c_double
 second=0.0_c_double
 do t=50,300
  growth=log(output(t))-log(output(t-1))
  mean=mean+growth
  second=second+growth*growth
 end do
 mean=mean/251.0_c_double
 second=second/251.0_c_double
 variance=second-mean*mean
end subroutine
end module
