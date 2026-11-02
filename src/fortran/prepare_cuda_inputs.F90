! Construct parameter-dependent CUDA inputs only; no policy solve or moments.
program prepare_cuda_inputs
use params
use solve, only: discretize, outputnetlabor
use, intrinsic :: iso_c_binding, only: c_double
use base_lib, only: boxmuller
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer :: gpu_jobs,gpu_io,gpu_unit,gpu_i,gpu_offset,ct,t,firmct,Act,Aprimect,Asavect,zct,zprimect,zsavect,kct,kprimect
character(len=16384) :: gpu_line
character(len=4096) :: gpu_manifest,gpu_dir,gpu_shared,gpu_path
character(len=64) :: gpu_extra
double precision :: gpu_x(GMMparam),runval,gpu_ac(knum,knum),gpu_base(knum,knum),gpu_ynet(Anum,znum,knum)
double precision, allocatable :: upifirmshocks1(:,:),upifirmshocks2(:,:)
logical :: gpu_server
double precision, volatile :: gpu_depreciated
double precision :: gpu_investment
interface
 function gpu_fma(a,b,c) bind(C,name="fma") result(value)
  import c_double
  real(c_double), value :: a,b,c
  real(c_double) :: value
 end function
end interface

if (command_argument_count()/=2) error stop 'Expected manifest and shared-output directory'
call get_command_argument(1,gpu_manifest)
call get_command_argument(2,gpu_shared)
gpu_server=trim(gpu_manifest)=='--serve'
allocate(pr_mat_z(znum,znum),pr_mat_z_theta(znum,znum,znum),pr_mat_A(Anum,Anum),&
 pr_mat_A_theta(Anum,Anum,Anum),pr_mat_exog(exognum,exognum),&
 pr_mat_exog_theta(exognum,exognum),A0(Anum),z0(znum),k0(knum),b0(bnum),&
 exogind(exognum,4),exogval(exognum,4),endogind(endognum,2),endogval(endognum,2),&
 ergdistA(Anum),spread0(spreadnum),Asimshocks(Tsim),Asimpos(Tsim),Asim(Tsim),&
 zfirmshocks(Tsim,Nfirm),zfirmpos(Tsim,Nfirm),pifirmshocks(Tsim,Nfirm),&
 upifirmshocks1(Tsim,Nfirm),upifirmshocks2(Tsim,Nfirm))
open(13,file="Asimshocks.txt")
do t=1,Tsim
	read(13,*) Asimshocks(t)
end do !t
close(13)

open(13,file="zfirmshocks.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  read(13,*) zfirmshocks(t,firmct)
  end do !firmct
end do !t
close(13)

open(13,file="upifirmshocks1.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  read(13,*) upifirmshocks1(t,firmct)
  end do !firmct
end do !t
close(13)

open(13,file="upifirmshocks2.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  read(13,*) upifirmshocks2(t,firmct)
  end do !firmct
end do !t
close(13)



do ct=1,Nfirm
 call boxmuller(upifirmshocks1(:,ct),upifirmshocks2(:,ct),pifirmshocks(:,ct),Tsim)
end do
open(newunit=gpu_unit,file=trim(gpu_shared)//'/RCC_pifirmshocks.bin',access='stream',&
 form='unformatted',status='replace',action='write')
write(gpu_unit) pifirmshocks
close(gpu_unit)
if (gpu_server) then
 write(*,'(A)') '{"family":"ready","protocol":1}'
 flush(6)
end if
do
 if (gpu_server) then
  read(*,'(A)',iostat=gpu_io) gpu_manifest
  if(gpu_io/=0 .or. trim(gpu_manifest)=='QUIT') exit
 end if
open(newunit=gpu_jobs,file=trim(gpu_manifest),status='old',action='read')
do
 read(gpu_jobs,'(A)',iostat=gpu_io) gpu_line
 if (gpu_io<0) exit
 if (gpu_io/=0) error stop 'Cannot read input manifest'
 read(gpu_line,*,iostat=gpu_io) gpu_dir,gpu_x
 if (gpu_io/=0 .or. .not.all(ieee_is_finite(gpu_x))) error stop 'Invalid parameter record'
#ifdef MODEL_DE
 theta=gpu_x(1)
 gpu_offset=1
#else
 gpu_offset=0
#endif
 rhoz=gpu_x(1+gpu_offset);sigmaz=gpu_x(2+gpu_offset);etak=gpu_x(3+gpu_offset)
 phi=gpu_x(4+gpu_offset);gamma=gpu_x(5+gpu_offset);sigmaA=gpu_x(6+gpu_offset)
 etaf=gpu_x(7+gpu_offset);etad=gpu_x(8+gpu_offset);sigmapi=gpu_x(9+gpu_offset)
 if (abs(rhoz)>=1.0d0 .or. sigmaz<=0.0d0 .or. sigmaA<=0.0d0) error stop 'Invalid process parameters'
 call discretize()
!now, need to compare the U(0,1) shocks to the A grid points
t = 1
Asimpos(t) = Anum/2
Asim(t) = A0(Asimpos(t))
do t=2,Tsim

	Act = Asimpos(t-1)
	Asavect=0
	runval = 0.0
	do Aprimect=1,Anum
		runval = runval+pr_mat_A(Act,Aprimect)
		if ((Asimshocks(t)<=runval).and.(Asavect==0)) then
			Asavect = Aprimect
		end if
	end do

	Asimpos(t) = Asavect
	Asim(t) = A0(Asavect)
end do !t

!now, do the exogenous firm-level simulation
t = 1
zfirmpos(t,:) = znum/2

do firmct=1,Nfirm
do t=2,Tsim

	zct = zfirmpos(t-1,firmct)
	zsavect=0
	runval = 0.0
	do zprimect=1,znum
		runval = runval+pr_mat_z(zct,zprimect)
		if ((zfirmshocks(t,firmct)<=runval).and.(zsavect==0)) then
			zsavect = zprimect
		end if
	end do

	zfirmpos(t,firmct) = zsavect
end do !t
end do !t


 do Act=1,Anum
  do zct=1,znum
   do kct=1,knum
    gpu_ynet(Act,zct,kct)=outputnetlabor(A0(Act),z0(zct),k0(kct))
   end do
  end do
 end do
 do kct=1,knum
  do kprimect=1,knum
   ! Preserve the screening cache's rounded investment subtraction, followed
   ! by its fused capital-cost combination. Exact action rescoring is separate.
   gpu_depreciated=(1.0-delta)*k0(kct)
   gpu_investment=k0(kprimect)-gpu_depreciated
   if(k0(kct)>0.0d0) then
    gpu_ac(kct,kprimect)=(etak/2.0)*(gpu_investment/k0(kct))**2.0*k0(kct)
   else
    gpu_ac(kct,kprimect)=0.0d0
   end if
   gpu_base(kct,kprimect)=gpu_fma(-(1.0-tau),gpu_ac(kct,kprimect),-invprice*k0(kprimect))
  end do
 end do
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_pr_mat_exog_theta.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) pr_mat_exog_theta
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_pr_mat_exog.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) pr_mat_exog
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_pr_mat_z.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) pr_mat_z
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_k0.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) k0
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_b0.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) b0
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_A0.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) A0
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_z0.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) z0
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_Asimpos.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) Asimpos
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_zfirmpos.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) zfirmpos
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_actionkbase.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) gpu_base
 close(gpu_unit)
 open(newunit=gpu_unit,file=trim(gpu_dir)//'/RCC_ynetcache.bin',access='stream',&
 form='unformatted',status='replace',action='write')
 write(gpu_unit) gpu_ynet
 close(gpu_unit)
end do
close(gpu_jobs)
 if (.not.gpu_server) exit
 write(*,'(A)') '{"family":"batch_complete","status":0}'
 flush(6)
end do
end program prepare_cuda_inputs
