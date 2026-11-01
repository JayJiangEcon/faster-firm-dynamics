! Command-line evaluator linked against either released or refactored RE solve.
! All economic calculations, including Ysim, remain inside the linked modules.
program evaluate_re
  use params
  use solve
  implicit none

  integer :: i, read_status
  character(len=128) :: arg
  double precision :: x(GMMparam), fx

  if (command_argument_count()/=GMMparam) then
    write(*, '(A,I0,A)') 'Expected ', GMMparam, ' RE parameters.'
    error stop 2
  end if

  do i=1,GMMparam
    call get_command_argument(i,arg)
    read(arg,*,iostat=read_status) x(i)
    if (read_status/=0) then
      write(*, '(A,I0,A,A)') 'Could not parse parameter ',i,': ',trim(arg)
      error stop 2
    end if
  end do

  fx = fGMM(x)
  write(*, '(A,ES24.16)') 'RCC_FGMM=', fx
  write(*, '(A,ES24.16)') 'RCC_RAW_Q=', GMMobjval
  write(*, '(A,17(ES24.16,1X))') 'RCC_MOMENTS=', GMMmodelmom
end program evaluate_re
