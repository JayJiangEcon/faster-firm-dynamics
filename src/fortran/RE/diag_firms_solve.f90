module params
! RCC_RE_STATE_COLLAPSED: theta=0 sufficient state is (A,z); lag columns are compatibility aliases.
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Firms w/Risky Debt + Agg Shocks
!
! Module to store parameter values
!
! Stephen Terry
!
! Original version: 11/7/2018
! + This Version: 2/7/2019 removes distributional IRFs, adds doVFI,
! switches to distributional simulation
! + This version: 2/14/2019 adds double IRF
! + This version: 4/23/2019 adds exact data matching
! + This version: 5/13/2019 adds issuance weighted spreads
! + This version: 4/14/2020 adds list of 25 micro + macro moments, k = 0 state
! + This version: 8/17/2020 allows for equity issuance, changes to capital only recovery, allows VD continuation
!!!!!!!!!!!!!!!!!!!!!!!!!!

!model parameters fixed outside any estimation
double precision, parameter :: alpha = 0.25 !capital revenue elasticity
double precision, parameter :: nu = 0.5 !labor revenue elasticity
double precision, parameter :: delta = 0.1 !capital depreciation rate
double precision, parameter :: invprice = 1.0 !relative price of investment
double precision, parameter :: rhoA = 0.95 !persistence of agg prod
double precision, parameter :: w = 0.5 !fixed wage
double precision, parameter :: r = 0.04 !fixed real risk free rate
double precision, parameter :: beta = 1.0/(1.0+r) !discount rate
double precision, parameter :: tau = 0.2 !fixed corporate income tax rate
double precision, parameter :: etab = 0.0 !quadratic debt expansion cost
!double precision, parameter :: etaf = 0.01 !fixed cost of equity issuance
!double precision, parameter :: etad = 0.03 !linear cost of equity issuance
double precision, parameter :: gammapi = 0.0 !MA(1) pi measurement error term
double precision, parameter :: theta = 0.0 !diagnosticity parameter

!model parameters that may vary with an estimation
double precision :: rhoz,sigmaz,etak,phi,gamma,sigmaA,sigmapi,etaf,etad

!GMM setup
integer, parameter :: GMMparam = 9 !number of parameters in estimation
integer, parameter :: GMMmom_small = 16 !number of moments in estimation
integer, parameter :: GMMmom = 17 !total number of moments in moment vectors and matrices

integer, parameter :: GMMmomselect(GMMmom_small) = (/ 2, 3, 4, 5, 6, 7, 8, 9, 10, &
  11, 12, 13, 14, 15, 16, 17 /)

integer, parameter :: numstats = 13 !number of statistics in moment calcs at the micro level
integer, parameter :: numstatsmacro = 3 !number of statistics in moment calcs at the micro level
integer, parameter :: verbose = 0 !write files if = 1, don't if = 0

double precision :: GMMdatawgt(GMMmom,GMMmom) !moment weighting matrix from the data
double precision :: GMMdatamom(GMMmom) !moment vector from the data

integer :: solveFLAG,bunchFLAG,distFLAG

integer, parameter :: knum = 30 !number of capital grid points for policies
integer, parameter :: bnum = 30 !number of debt grid points
integer, parameter :: znum = 5 !number of idio prod grid points
integer, parameter :: anum = 5 !number of agg prod grid points
integer, parameter :: exognum= anum*znum !RE sufficient exog state grid points (A,z)
integer, parameter :: endognum = knum*bnum !number of total endog state grid points (k,b)

double precision, parameter :: nsigmaz = 2.5 !multiple of z std dev's to cover
double precision, parameter :: nsigmaA = 1.25 !multiple of A std dev's to cover
double precision, parameter :: kmin = 0.75 !min k level
double precision, parameter :: kmax = 7.5 !max k level
double precision, parameter :: bmin = 0.75 !min b level
double precision, parameter :: bmax = 7.5 !max b level

integer, parameter :: spreadnum = 250 !number of grid points for spread tracking grid
double precision, parameter :: spreadmin = 0.0 !min of spread tracking grid
double precision, parameter :: spreadmax = 1.0 !max of spread tracking grid

integer :: dfmaxit = 50 !max number of default iterations for convergence
integer :: polmaxit = 20 !max number of policy iterations for convergence
integer :: accelmaxit = 100 !number of pushforwards of the value function
integer :: distmaxit = 1500 !max number of dist pushforwards
integer, parameter :: solwindow = 10 !number of grid points in either direction for policies
double precision :: disttol = 1.0e-6 !convergence tolerance on distribution
double precision, parameter :: vftol = 1.0e-4 !convergence tolerance on vf
double precision :: poltol = 1.0e-4 !convergence tolerance on default policies

integer, parameter :: doVFI = 1 !if 1, do policy iteration, if 0 read in the policy and value files

integer, parameter :: doRATprice = 0 !if 1, do rational pricing, if 0, pricing matches theta parameter
integer, parameter :: doRATfirms = 0 !if 1, do rational firm expectations, if 0 firm expectations matches theta parameter

double precision, parameter :: junkspread = 0.025 !spread level which is defined as junk

!simulation parameters
integer, parameter :: seedint = 25477 !seed for random shocks
integer, parameter :: Terg = 50 !number of periods to burn in
integer, parameter :: Tgood = 250 !number of periods to use for actual simulation & fcst updates
integer, parameter :: Tsim = Terg + Tgood + 2 !total number of periods to simulate
integer, parameter :: Trepresent = 10 !period in which to extract representative dist
integer, parameter :: Nfirm = 1000 !number of firms to simulate
integer, parameter :: Tfirm = 6 !periods over which to compute moments at the firm level
integer, parameter :: doREADdraws = 1 !if 1, read draws from existing files, if 0 generate them


!IRF parameters
integer, parameter :: doIRF = 0 !if 1 do IRF, if 0 don't do IRF
integer, parameter :: Tirf = 75 !total number of periods to simulate for each economy
integer, parameter :: Tstartirf = 50 !period in which shock sequence starts
integer, parameter :: Tnumirf = 2 !number of periods with the shock sequence
integer, parameter :: seqAposirf(Tnumirf) = (/ 1,-1 /) !the path of agg prod to feed into each economy
integer, parameter :: Nirf = 10000 !number of economies to simulate
!integer, parameter :: Nirf = 200 !number of economies to simulate

double precision, parameter :: sizeIRFshock(Tnumirf) = (/1.5*0.006377651836837,-1.0*0.006377651836837/) !the average size of the shocks
double precision :: probIRFshock(Tnumirf)
double precision, parameter :: IRFtol = 1.0e-4
integer, parameter :: IRFmaxit = 100

!Data matching parameters
integer, parameter :: doMATCH = 0 !if 1, do data match, if 0 don't do data match
integer, parameter :: doreadMATCH = 0 !if 1, read data match
integer, parameter :: Tmatchser = 12 !length of the series of growth rates you desire to match
integer, parameter :: Tmatch = Tmatchser + 1 !total length of series in data matching
integer, parameter :: Tmatchsimstart = 1441 !period in the unconditional simulation to initialize distribution with
double precision, parameter :: dIdata(Tmatch) = (/ &
0.0,&
0.0389,&
-0.0207,&
-0.1680,&
0.0153,&
0.0529,&
0.0579,&
0.0044,&
0.0295,&
-0.0200,&
-0.0331,&
0.0125,&
0.0282 /)

integer, parameter :: maxAiters = 150 !number of iterations on clearing level of A
double precision, parameter :: Atol = 1.0e-4 !tolerance on convergence of A guess

!other model quantities which need to be globally available
double precision :: GMMobjval,GMMmodelmom(GMMmom)

integer :: seeddim,GEct

double precision, allocatable :: Vtheta(:,:),Vthetaold(:,:),Vndtheta(:,:),pr_mat_z(:,:),&
	pr_mat_z_theta(:,:,:),pr_mat_A(:,:),pr_mat_A_theta(:,:,:),&
	pr_mat_exog(:,:),A0(:),z0(:),kpol(:,:),kpolold(:,:),exogval(:,:),k0(:),b0(:),&
	pr_mat_exog_theta(:,:),V(:,:),Vold(:,:),Vnd(:,:),bpol(:,:),bpolold(:,:),endogval(:,:),&
	dfpol(:,:),qmat(:,:),dist(:,:),distold(:,:),distk(:),distb(:),&
	Ysim(:),Ksim(:),Isim(:),qsim(:),Bsim(:),Asimshocks(:),Asim(:),dfsim(:),distsim(:,:,:),&
  BtoKsim(:),YsimIRF(:,:,:),KsimIRF(:,:,:),IsimIRF(:,:,:),qsimIRF(:,:,:),&
  BsimIRF(:,:,:),AsimshocksIRF(:,:),AsimIRF(:,:,:),dfsimIRF(:,:,:),&
  BtoKsimIRF(:,:,:),distbeforeIRF(:,:),distafterIRF(:,:),irateheterosim(:,:),&
  Dbheterosim(:,:),zfirmshocks(:,:),pimat(:,:),pifirm(:,:),pifcstfirm(:,:),&
  ifirm(:,:),bfirm(:,:),kfirm(:,:),fefirm(:,:),qfirm(:,:),Amatch(:),pAmatch(:),&
  Ymatch(:),Imatch(:),dImatch(:),Kmatch(:),qmatch(:),dfmatch(:),Bmatch(:),&
  BtoKmatch(:),distmatch(:,:,:),spreadsim(:),spreadsimIRF(:,:,:),&
  spreadmatch(:),dspreadmatch(:),pifirmshocks(:,:),retbondsim(:),spread0(:),&
  spreadtracksim(:,:),spreadBAAsim(:),spreadBAAmatch(:),spreadtrackmatch(:,:),&
  dspreadBAAmatch(:),Nsim(:),sdMPKsim(:),sdMPKmatch(:),Nmatch(:),dfYsim(:),&
  dfBsim(:),dfIsim(:),issjunksim(:),bspreadtracksim(:,:),ergdistA(:),imposeshocksIRF(:,:),&
  divposmat(:,:),Vdtheta(:),Vdthetaold(:),FEsim(:),FEsimIRF(:,:,:),Pisim(:),PisimIRF(:,:,:),&
  Pifsim(:),PifsimIRF(:,:,:),sprfirm(:,:)

integer, allocatable :: kpolind(:,:),kpolindold(:,:),exogind(:,:),bpolind(:,:),&
	bpolindold(:,:),endogind(:,:),kbpolind(:,:),kbpolindold(:,:),dfpolind(:,:),dfpolindold(:,:),&
	Asimpos(:),seedarray(:),AsimposIRF(:,:,:),zfirmpos(:,:),endogfirmpos(:,:),&
  dfirmpos(:,:),kfirmpos(:,:),bfirmpos(:,:),Amatchpos(:)

end module params

module solve
use params
use base_lib
use omp_lib
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Firms w/Risky Debt + Agg Shocks
!
! Module to store functions and subroutines for solution of the model
!
! Stephen Terry
!
! This version: 11/7/2018
!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine configure_fidelity()
use params
implicit none
character(len=64) :: value
integer :: status

dfmaxit = 50
polmaxit = 20
accelmaxit = 100
distmaxit = 1500
disttol = 1.0d-6
poltol = 1.0d-4

call get_environment_variable("RCC_DFMAXIT",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) dfmaxit
call get_environment_variable("RCC_POLMAXIT",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) polmaxit
call get_environment_variable("RCC_ACCELMAXIT",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) accelmaxit
call get_environment_variable("RCC_DISTMAXIT",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) distmaxit
call get_environment_variable("RCC_DISTTOL",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) disttol
call get_environment_variable("RCC_POLTOL",value,status=status)
if (status==0 .and. len_trim(value)>0) read(value,*) poltol

if (dfmaxit<2 .or. polmaxit<1 .or. accelmaxit<1 .or. distmaxit<1 &
    .or. disttol<=0.0d0 .or. poltol<=0.0d0) then
  error stop "invalid RCC fidelity controls"
end if
write(*,'(A,4(I0,1X),2(ES12.4,1X))') "RCC_FIDELITY_CONTROLS=", &
  dfmaxit,polmaxit,accelmaxit,distmaxit,disttol,poltol
end subroutine configure_fidelity


double precision function fGMM(x)
use params
use base_lib
use omp_lib
implicit none

!declarations
integer :: ct
double precision :: x(GMMparam),start,werr,stage_start
double precision :: sigma_profile_m0(GMMmom),sigma_profile_mp(GMMmom),&
  sigma_profile_mm(GMMmom),sigma_profile_c1(GMMmom),sigma_profile_c2(GMMmom),&
  sigma_profile_q0
integer :: sigma_profile_status
character(len=32) :: sigma_profile_mode

call configure_fidelity()

! OMP_NUM_THREADS sets the thread count for this evaluation.

!$omp parallel
write(*,*) "Parallel hello!"
!$omp end parallel

start = omp_get_wtime()
write(*,*) "#############################################"
write(*,*) "#############################################"
write(*,*) "Starting model solution"
write(*,*) "#############################################"
write(*,*) "#############################################"
write(*,*) " "

!read in the parameters
rhoz = x(1)
sigmaz = x(2)
etak = x(3)
phi = x(4)
gamma = x(5)
sigmaA = x(6)
etaf = x(7)
etad = x(8)
sigmapi = x(9)

!write some parameters
write(*,*) "##### Parameters for this solution "
write(*,*) "rhoz = ",rhoz
write(*,*) "sigmaz = ",sigmaz
write(*,*) "etak = ",etak
write(*,*) "phi = ",phi
write(*,*) "gamma = ",gamma
write(*,*) "sigmaA = ",sigmaA
write(*,*) "etaf = ",etaf
write(*,*) "etad = ",etad
write(*,*) "sigmapi = ",sigmapi
write(*,*) " "

!allocate some of the globally available arrays
allocate(Vtheta(exognum,endognum),Vthetaold(exognum,endognum),pr_mat_z(znum,znum),&
	pr_mat_z_theta(znum,znum,znum),pr_mat_A(Anum,Anum),pr_mat_A_theta(Anum,Anum,Anum),&
	pr_mat_exog(exognum,exognum),A0(Anum),z0(znum),kpol(exognum,endognum),&
	kpolold(exognum,endognum),exogval(exognum,4),k0(knum),b0(bnum),&
	pr_mat_exog_theta(exognum,exognum),V(exognum,endognum),Vold(exognum,endognum),&
	bpol(exognum,endognum),bpolold(exognum,endognum),endogval(endognum,2),&
	kpolind(exognum,endognum),kpolindold(exognum,endognum),exogind(exognum,4),&
	bpolind(exognum,endognum),bpolindold(exognum,endognum),endogind(endognum,2),&
	kbpolind(exognum,endognum),kbpolindold(exognum,endognum),Vndtheta(exognum,endognum),&
	Vnd(exognum,endognum),dfpol(exognum,endognum),dfpolind(exognum,endognum),&
	dfpolindold(exognum,endognum),qmat(exognum,endognum),&
	dist(exognum,endognum),distold(exognum,endognum),distk(knum),distb(bnum),&
  Ysim(Tsim),Ksim(Tsim),Isim(Tsim),qsim(Tsim),Bsim(Tsim),Asimpos(Tsim),Asimshocks(Tsim),&
  Asim(Tsim),dfsim(Tsim),distsim(exognum,endognum,Tsim),BtoKsim(Tsim),&
  YsimIRF(Tirf,Nirf,2),KsimIRF(Tirf,Nirf,2),IsimIRF(Tirf,Nirf,2),qsimIRF(Tirf,Nirf,2),&
  BsimIRF(Tirf,Nirf,2),AsimshocksIRF(Tirf,Nirf),AsimIRF(Tirf,Nirf,2),&
  dfsimIRF(Tirf,Nirf,2),BtoKsimIRF(Tirf,Nirf,2),AsimposIRF(Tirf,Nirf,2),&
  distbeforeIRF(exognum,endognum),distafterIRF(exognum,endognum),&
  irateheterosim(Tsim,6),Dbheterosim(Tsim,6),zfirmshocks(Tsim,Nfirm),&
  zfirmpos(Tsim,Nfirm),endogfirmpos(Tsim,Nfirm),dfirmpos(Tsim,Nfirm),&
  pimat(exognum,endognum),kfirmpos(Tsim,Nfirm),bfirmpos(Tsim,Nfirm),&
  pifirm(Tsim,Nfirm),pifcstfirm(Tsim,Nfirm),ifirm(Tsim,Nfirm),&
  bfirm(Tsim,Nfirm),kfirm(Tsim,Nfirm),fefirm(Tsim,Nfirm),qfirm(Tsim,Nfirm),&
  Amatch(Tmatch),pAmatch(Tmatch),Ymatch(Tmatch),Imatch(Tmatch),dImatch(Tmatch),&
  Kmatch(Tmatch),qmatch(Tmatch),dfmatch(Tmatch),Bmatch(Tmatch),&
  BtoKmatch(Tmatch),Amatchpos(Tmatch),distmatch(exognum,endognum,Tmatch),&
  spreadsim(Tsim),spreadsimIRF(Tirf,Nirf,2),spreadmatch(Tmatch),&
  dspreadmatch(Tmatch),pifirmshocks(Tsim,Nfirm),retbondsim(Tsim),spread0(spreadnum),&
  spreadtracksim(spreadnum,Tsim),spreadBAAsim(Tsim),spreadBAAmatch(Tmatch),&
  spreadtrackmatch(spreadnum,Tmatch),dspreadBAAmatch(Tmatch),Nsim(Tsim),&
  sdMPKsim(Tsim),sdMPKmatch(Tmatch),Nmatch(Tmatch),dfYsim(Tsim),dfBsim(Tsim),&
  dfIsim(Tsim),issjunksim(Tsim),bspreadtracksim(spreadnum,Tsim),ergdistA(Anum),&
  imposeshocksIRF(Tirf,Nirf),divposmat(exognum,endognum),Vdtheta(exognum),Vdthetaold(exognum),&
  FEsim(Tsim),FEsimIRF(Tirf,Nirf,2),Pisim(Tsim),PisimIRF(Tirf,Nirf,2),Pifsim(Tsim),&
    PifsimIRF(Tirf,Nirf,2),sprfirm(Tsim,Nfirm))

!read in GMM moment weight matrix & moment vector
GMMdatawgt(:,:) = 0.0
open(13,file="_GMMdatawgt.txt")
do ct=1,GMMmom
read(13,*) GMMdatawgt(ct,:)
end do !ct
close(13)

GMMdatamom(:) = 0.0
open(13,file="_GMMdatamom.txt")
do ct=1,GMMmom
read(13,*) GMMdatamom(ct)
end do !ct
close(13)

!set up the model grids and discretize the driving process
stage_start = omp_get_wtime()
call discretize()
write(*, '(A,ES24.16)') 'RCC_TIME_DISCRETIZE=',omp_get_wtime()-stage_start

!now, do the policy iteration
stage_start = omp_get_wtime()
call politeration()
write(*, '(A,ES24.16)') 'RCC_TIME_POLICY=',omp_get_wtime()-stage_start

!Policy iteration already leaves the distribution at the final policy.  Only
!the read-policy path still needs a distribution call here.
stage_start = omp_get_wtime()
if (doVFI==0) call distribution()
write(*, '(A,ES24.16)') 'RCC_TIME_DISTRIBUTION=',omp_get_wtime()-stage_start

!now, do the simulation
stage_start = omp_get_wtime()
call get_environment_variable('RCC_SIGMAPI_PROFILE',sigma_profile_mode,&
  status=sigma_profile_status)
if (sigma_profile_status==0 .and. trim(sigma_profile_mode)=='1') then
  if (sigmapi/=0.0d0) error stop 'RE sigma profile requires x(9)=0'
  call simulate(.true.,sigma_profile_m0,sigma_profile_mp,sigma_profile_mm)
  sigma_profile_q0 = GMMobjval
  sigma_profile_c1 = 0.5d0*(sigma_profile_mp-sigma_profile_mm)
  sigma_profile_c2 = 0.5d0*(sigma_profile_mp+sigma_profile_mm)-sigma_profile_m0
  write(*,'(A,17(ES24.16,1X))') 'RCC_SIGMAPI_M0=',sigma_profile_m0
  write(*,'(A,17(ES24.16,1X))') 'RCC_SIGMAPI_MP=',sigma_profile_mp
  write(*,'(A,17(ES24.16,1X))') 'RCC_SIGMAPI_MM=',sigma_profile_mm
  write(*,'(A,17(ES24.16,1X))') 'RCC_SIGMAPI_C1=',sigma_profile_c1
  write(*,'(A,17(ES24.16,1X))') 'RCC_SIGMAPI_C2=',sigma_profile_c2
  GMMmodelmom = sigma_profile_m0
  GMMobjval = sigma_profile_q0
else
  call simulate(.false.,sigma_profile_m0,sigma_profile_mp,sigma_profile_mm)
end if
write(*, '(A,ES24.16)') 'RCC_TIME_SIMULATION=',omp_get_wtime()-stage_start

!IRF simulations
if (doIRF==1) then
call simulateIRF()
end if

!write some output files
if (verbose==1) then
  call writefiles()
end if

fGMM = GMMobjval

open(13,file="GMMobj.txt",position="append")
write(13,*) fGMM
close(13)

open(13,file="GMMparam.txt",position="append")
write(13,*) x
close(13)

open(13,file="GMMmodelmom.txt",position="append")
write(13,*) GMMmodelmom
close(13)


open(13,file="GMMdiag.txt",position="append")
write(13,*) omp_get_wtime()-start,solveFLAG,bunchFLAG,distFLAG
close(13)

!deallocate to avoid memory clutter
deallocate(Vtheta,Vthetaold,pr_mat_z,pr_mat_z_theta,pr_mat_A,pr_mat_A_theta,&
	pr_mat_exog,A0,z0,kpol,kpolold,exogval,k0,b0,pr_mat_exog_theta,V,Vold,&
	bpol,bpolold,endogval,kpolind,kpolindold,exogind,bpolind,bpolindold,endogind,&
	kbpolind,kbpolindold,Vndtheta,Vnd,dfpol,dfpolind,dfpolindold,qmat,&
	dist,distold,distk,distb,Ysim,Ksim,Isim,qsim,Bsim,Asimpos,Asimshocks,Asim,&
  dfsim,distsim,BtoKsim,YsimIRF,KsimIRF,IsimIRF,qsimIRF,BsimIRF,AsimshocksIRF,&
  AsimIRF,dfsimIRF,BtoKsimIRF,AsimposIRF,distbeforeIRF,distafterIRF,&
  irateheterosim,Dbheterosim,zfirmshocks,zfirmpos,endogfirmpos,&
  dfirmpos,pimat,kfirmpos,bfirmpos,pifirm,pifcstfirm,ifirm,bfirm,kfirm,fefirm,&
  qfirm,Amatch,pAmatch,Ymatch,Imatch,dImatch,Kmatch,qmatch,dfmatch,Bmatch,&
  BtoKmatch,Amatchpos,distmatch,spreadsim,spreadsimIRF,spreadmatch,&
  dspreadmatch,pifirmshocks,retbondsim,spread0,spreadtracksim,spreadBAAsim,&
  spreadBAAmatch,spreadtrackmatch,dspreadBAAmatch,Nsim,sdMPKsim,&
  sdMPKmatch,Nmatch,dfYsim,dfBsim,dfIsim,issjunksim,bspreadtracksim,ergdistA,&
  imposeshocksIRF,divposmat,Vdtheta,Vdthetaold,FEsim,FEsimIRF,Pisim,PisimIRF,&
  Pifsim,PifsimIRF,sprfirm)

write(*,*) "Done with function evaluation at ",omp_get_wtime()-start," secs."
write(*,*) "#############################################"

end function fGMM


subroutine simulateIRF()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,t,Asavect,Act,Aprimect,firmct,zsavect,zct,zprimect,&
	endogct,exogct,Amin1ct,zmin1ct,dfct,kprimect,bprimect,kbpolct,&
	endogprimect,exogprimect,simct,shockct,irfit

double precision :: runval,dumkvec(knum),dumbvec(bnum),yval,zval,Aval,&
	kval,ival,bval,qval,kprimeval,bprimeval,issdenomval,poskmass,poskprimemass,&
  chgvec(Nirf),pa,pb,realizedval,desiredval,irferr

write(*,*) "##### Doing model IRF simulation "
write(*,*) " "

!this block seeds the random draws
call random_seed(size=seeddim)
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint+ct
end do !ct
call random_seed(put=seedarray)

!call the random draws for the productivity process
call random_number(AsimshocksIRF); !U(0,1) in each cell
call random_number(imposeshocksIRF); !U(0,1) in each cell

deallocate(seedarray)

!!!!!!!!!!
!!!!do the purely exogenous agg prod sim
!!!!!!!!!!

!now, need to compare the U(0,1) shocks to the A grid points
do simct=1,Nirf

  t = 1
  AsimposIRF(t,simct,:) = Anum/2
  AsimIRF(t,simct,:) = A0(AsimposIRF(t,simct,:))

  !start off with normal unconditional simulation (for both)
  do t=2,Tstartirf-1

  	Act = AsimposIRF(t-1,simct,1)
  	Asavect=0
  	runval = 0.0
  	do Aprimect=1,Anum
  		runval = runval+pr_mat_A(Act,Aprimect)
  		if ((AsimshocksIRF(t,simct)<=runval).and.(Asavect==0)) then
  			Asavect = Aprimect
  		end if
  	end do

  	AsimposIRF(t,simct,:) = Asavect
  	AsimIRF(t,simct,:) = A0(Asavect)
  end do !t

  !continue with normal unconditional simulation (for 2)
  do t=Tstartirf,Tirf

    Act = AsimposIRF(t-1,simct,2)
    Asavect=0
    runval = 0.0
    do Aprimect=1,Anum
      runval = runval+pr_mat_A(Act,Aprimect)
      if ((AsimshocksIRF(t,simct)<=runval).and.(Asavect==0)) then
        Asavect = Aprimect
      end if
    end do

    AsimposIRF(t,simct,2) = Asavect
    AsimIRF(t,simct,2) = A0(Asavect)
  end do !t

end do !simct

  !then feed in the desired sequence of shocks
  do t=1,Tnumirf
    pa = 0.0
    pb = 1.0


    do irfit = 1,IRFmaxit


          probIRFshock(t) = 0.5*(pa+pb)

      do simct=1,Nirf

          !impose shock?
          if (imposeshocksIRF(Tstartirf+t-1,simct)<probIRFshock(t)) then


            AsimposIRF(Tstartirf+t-1,simct,1) = AsimposIRF(Tstartirf+t-2,simct,1)+seqAposirf(t)

          !or draw normally?
          else

            Act = AsimposIRF(Tstartirf+t-2,simct,1)
            Asavect=0
            runval = 0.0
            do Aprimect=1,Anum
              runval = runval+pr_mat_A(Act,Aprimect)
              if ((AsimshocksIRF(t,simct)<=runval).and.(Asavect==0)) then
                Asavect = Aprimect
              end if
            end do

            AsimposIRF(Tstartirf+t-1,simct,1) = Asavect

          end if

          if (AsimposIRF(Tstartirf+t-1,simct,1)<1) AsimposIRF(Tstartirf+t-1,simct,1) = 1
          if (AsimposIRF(Tstartirf+t-1,simct,1)>=Anum) AsimposIRF(Tstartirf+t-1,simct,1) = Anum

          AsimIRF(Tstartirf+t-1,simct,1) = A0(AsimposIRF(Tstartirf+t-1,simct,1))
      end do !simct

      chgvec(:) = 0.0
      do simct=1,Nirf

        chgvec(simct) = log(A0(AsimposIRF(Tstartirf+t-1,simct,1)))-log(A0(AsimposIRF(Tstartirf+t-1,simct,2)))
      end do !simct


  write(*,*) "guessed IRF shock prob ",probIRFshock(t)

  !what is realized diff between shock and no shock economies?
  realizedval = sum(chgvec)/dble(Nirf)
  write(*,*) "realized val",realizedval

  !what is actual diff?
  if (t==1) then
    desiredval = sizeIRFshock(1)
  else if (t==2) then
    desiredval = sum(sizeIRFshock)
  end if
  write(*,*) "desired val",desiredval

  !what is error?
  irferr = realizedval - desiredval

  write(*,*) "irf err = ",irferr
  write(*,*) " "

  if (abs(irferr)<IRFtol) then
    exit
  else
    if (seqAposirf(t)>0) then

      if (irferr>0.0) then
        pb = probIRFshock(t)
      else
        pa = probIRFshock(t)
      end if
    else
      if (irferr<0.0) then
        pb = probIRFshock(t)
      else
        pa = probIRFshock(t)
      end if
    end if
  end if

  end do !irf it


end do !t

do simct=1,Nirf

  !then, simulate the remaining path unconditionally
  do t=Tstartirf+Tnumirf,Tirf

  	Act = AsimposIRF(t-1,simct,1)
  	Asavect=0
  	runval = 0.0
  	do Aprimect=1,Anum
  		runval = runval+pr_mat_A(Act,Aprimect)
  		if ((AsimshocksIRF(t,simct)<=runval).and.(Asavect==0)) then
  			Asavect = Aprimect
  		end if
  	end do

  	AsimposIRF(t,simct,1) = Asavect
  	AsimIRF(t,simct,1) = A0(Asavect)
  end do !t

end do !simct

!!!!!!!!!!
!!!!do the endogenous micro-level distributional simulation
!!!!!!!!!!
distbeforeIRF(:,:) = 0.0
distafterIRF(:,:) = 0.0

!now, initialize the agg variables
YsimIRF(:,:,:) = 0.0
KsimIRF(:,:,:) = 0.0
IsimIRF(:,:,:) = 0.0
qsimIRF(:,:,:) = 0.0
BsimIRF(:,:,:) = 0.0
dfsimIRF(:,:,:) = 0.0
BtoKsimIRF(:,:,:) = 0.0
spreadsimIRF(:,:,:) = 0.0
PisimIRF(:,:,:) = 0.0
PifsimIRF(:,:,:) = 0.0
FEsimIRF(:,:,:) = 0.0

!start with a representative period distribution
do simct=1,Nirf


  !with and without shock ct
  do shockct = 1,2

    !do the distributional initialization
    distbeforeIRF(:,:) = 1.0
    distbeforeIRF(:,:) = distbeforeIRF(:,:) /sum(distbeforeIRF(:,:))

  !now do the endogenous portions for subsequent periods
  do t=1,Tirf-1

    distafterIRF(:,:) = 0.0

    if (t==1) then
      write(*,*) "Doing period t = ",t," sim ",simct," shock = ",shockct
    end if

    if (t>1) then
      Amin1ct=AsimposIRF(t-1,simct,shockct)
    else if (t==1) then
      Amin1ct = Anum/2
    end if

    Act = AsimposIRF(t,simct,shockct)
    Aval = AsimIRF(t,simct,shockct)

    Aprimect=AsimposIRF(t+1,simct,shockct)


  do zmin1ct=1,1
  do zct=1,znum

    exogct = znum*(Act-1) + zct

    !loop over endog states
    do endogct=1,endognum
        if (distbeforeIRF(exogct,endogct)>0.0) then

        !does the firm default?
        dfct = dfpolind(exogct,endogct)

        if (dfct==0) then

    		  kbpolct = kbpolind(exogct,endogct)

      	else if (dfct==1) then

          kbpolct = 1

      	end if

        do zprimect=1,znum

          exogprimect = Anum*znum*znum*(Act-1)+znum*znum*(Aprimect-1) &
            + znum*(zct-1) + zprimect

          distafterIRF(exogprimect,kbpolct) = &
            distafterIRF(exogprimect,kbpolct) +&
            distbeforeIRF(exogct,endogct)*pr_mat_z(zct,zprimect)

        end do !zprimect

      end if
    end do !endogct

  end do !zct
  end do !zmin1ct

issdenomval = 0.0

  poskmass = 0.0
  poskprimemass = 0.0
  !now, actually do the insertion of the information into the aggregates
    do zct=1,znum
      do zmin1ct=1,1

        exogct = znum*(Act-1) + zct
        zval = z0(zct)

        do endogct=1,endognum
          !distribution >0?
          if (distbeforeIRF(exogct,endogct)>0.0) then

          kval = endogval(endogct,1)
          bval = endogval(endogct,2)

          dfct = dfpolind(exogct,endogct)

          if (dfct==0) then

            kbpolct = kbpolind(exogct,endogct)

          else if (dfct==1) then

            kbpolct = 1

          end if

          kprimeval = endogval(kbpolct,1)
          bprimeval = endogval(kbpolct,2)
        	yval = output(Aval,zval,kval)* (1.0-dble(dfct)*(1.0-gamma))
        	ival = kprimeval - (1.0-delta)*kval
        	!qval = qmat(exogct,endogct)
          qval = qmat(exogct,kbpolct)

          if (dfct==1) then
            yval = 0.0
          end if

          YsimIRF(t,simct,shockct) = YsimIRF(t,simct,shockct) + &
            yval*distbeforeIRF(exogct,endogct)
          KsimIRF(t,simct,shockct) = KsimIRF(t,simct,shockct) + &
            kval*distbeforeIRF(exogct,endogct)
        	IsimIRF(t,simct,shockct) = IsimIRF(t,simct,shockct) + &
            ival*distbeforeIRF(exogct,endogct)
        	qsimIRF(t,simct,shockct) = qsimIRF(t,simct,shockct) + &
            qval*distbeforeIRF(exogct,endogct)
        	BsimIRF(t,simct,shockct) = BsimIRF(t,simct,shockct) + &
            bval*distbeforeIRF(exogct,endogct)
        	dfsimIRF(t,simct,shockct) = dfsimIRF(t,simct,shockct) + &
            dble(dfct)*distbeforeIRF(exogct,endogct)
          if (kval>0.0) then
            BtoKsimIRF(t,simct,shockct) = BtoKsimIRF(t,simct,shockct) +&
              (bval/kval)*distbeforeIRF(exogct,endogct)
          end if
          spreadsimIRF(t,simct,shockct) = spreadsimIRF(t,simct,shockct) + &
            (1.0/qval-(1.0+r))*distbeforeIRF(exogct,endogct)
          issdenomval = issdenomval + bprimeval*kprimeval*&
            distbeforeIRF(exogct,endogct)

            if (kval>0.0) then
            PisimIRF(t,simct,shockct) = PisimIRF(t,simct,shockct) + &
              (pimat(exogct,endogct))*distbeforeIRF(exogct,endogct)
            poskmass = poskmass + distbeforeIRF(exogct,endogct)
            end if

            if (kprimeval>0.0) then
              do exogprimect=1,exognum
                PifsimIRF(t+1,simct,shockct) = PifsimIRF(t+1,simct,shockct) + (pimat(exogprimect,kbpolct))*&
                  pr_mat_exog_theta(exogct,exogprimect)*distbeforeIRF(exogct,endogct)
              end do !exogprimect
              poskprimemass = poskprimemass + distbeforeIRF(exogct,endogct)
            end if



          !end distribution selection point
          end if
        end do !endogct

      end do !zmin1ct
    end do !zct
    !PisimIRF(t,simct,shockct) = PisimIRF(t,simct,shockct)/poskmass
    !PifsimIRF(t+1,simct,shockct) = PifsimIRF(t+1,simct,shockct)/poskprimemass

    !round off and transfer the distributions to the next period
    distafterIRF(:,:) = distafterIRF(:,:)&
      /sum(distafterIRF(:,:))

    distbeforeIRF = distafterIRF

    !spreadsimIRF(t,simct,shockct) = spreadsimIRF(t,simct,shockct)/issdenomval


  end do !t
  end do !shockct

end do !simct
FEsimIRF = PisimIRF - PifsimIRF

write(*,*) " "

end subroutine simulateIRF



subroutine simulate(profile_sigma_mode,profile_m0,profile_mp,profile_mm)
use params
use base_lib
use omp_lib
implicit none

logical, intent(in) :: profile_sigma_mode
double precision, intent(out) :: profile_m0(GMMmom),profile_mp(GMMmom),profile_mm(GMMmom)

integer :: ct,t,Asavect,Act,Aprimect,firmct,zsavect,zct,zprimect,&
	endogct,exogct,Amin1ct,zmin1ct,dfct,kprimect,bprimect,kbpolct,&
	endogprimect,exogprimect,ct1,ct2,spreadind,found40,found60,spr0ct,&
  sprmedct,Nt(Tfirm),Ni(Nfirm),okobsflag,oldexogct,oldAct,oldAmin1ct

double precision :: runval,dumkvec(knum),dumbvec(bnum),yval,zval,Aval,&
	kval,ival,bval,qval,kprimeval,bprimeval,heterowgtvec(6),&
  dumvec(numstats),meanvec(numstats),covmat(numstats,numstats),&
    corrmat(numstats,numstats),pival,Amin1val,issdenomval,&
    upifirmshocks1(Tsim,Nfirm),upifirmshocks2(Tsim,Nfirm),&
    spreadval,spreadwgt,fracissuing,sum40,sum60,EMPK,EMPK2,nval,&
    junkdenomval,dumvecmacro(numstatsmacro),meanvecmacro(numstatsmacro),&
    covmatmacro(numstatsmacro,numstatsmacro),corrmatmacro(numstatsmacro,numstatsmacro),&
    momdiffvec(GMMmom_small),quadformvec(GMMmom_small),weightmat(GMMmom_small,GMMmom_small),&
    spreadvec(Nfirm),poskmass,poskprimemass,sprval,meant(Tfirm,numstats),meani(Nfirm,numstats)

write(*,*) "##### Doing model unconditional simulation "
write(*,*) " "

!this block seeds the random draws
if (doREADdraws==0) then
call random_seed(size=seeddim)
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint+ct
end do !ct
call random_seed(put=seedarray)

!call the random draws for the productivity process
call random_number(Asimshocks); !U(0,1) in each cell
call random_number(zfirmshocks); !U(0,1) in each cell
call random_number(upifirmshocks1); !U(0,1) in each cell
call random_number(upifirmshocks2); !U(0,1) in each cell
deallocate(seedarray)

open(13,file="Asimshocks.txt")
do t=1,Tsim
	write(13,*) Asimshocks(t)
end do !t
close(13)

open(13,file="zfirmshocks.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  write(13,*) zfirmshocks(t,firmct)
  end do !firmct
end do !t
close(13)

open(13,file="upifirmshocks1.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  write(13,*) upifirmshocks1(t,firmct)
  end do !firmct
end do !t
close(13)

open(13,file="upifirmshocks2.txt")
do t=1,Tsim
  do firmct=1,Nfirm
	  write(13,*) upifirmshocks2(t,firmct)
  end do !firmct
end do !t
close(13)

else if (doREADdraws==1) then


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


end if !doREADdraws

!convert uniform draws to normal draws for profit meas err

do ct=1,Nfirm
  call boxmuller(upifirmshocks1(:,ct),upifirmshocks2(:,ct),pifirmshocks(:,ct),Tsim)

end do !ct

!!!!!!!!!!
!!!!do the purely exogenous agg prod sim
!!!!!!!!!!

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

!!!!!!!!!!
!!!!do the endogenous micro-level distributional simulation
!!!!!!!!!!
distsim(:,:,:) = 0.0
pifirm(:,:) = 0.0
pifcstfirm(:,:) = 0.0
kfirm(:,:) = 0.0
ifirm(:,:) = 0.0
bfirm(:,:) = 0.0
qfirm(:,:) = 0.0
sprfirm(:,:) = 0.0
endogfirmpos(:,:) = 0
!zfirmpos(:,:) = 0
dfirmpos(:,:) = 0
kfirmpos(:,:) = 0
bfirmpos(:,:) = 0

!start with a flat distribution around current agg states
t = 1
Amin1ct = Anum/2
oldAmin1ct=Anum/2
Act = Asimpos(1)
oldAct = Asimpos(1)
do zmin1ct=1,1
  do zct=1,znum

    exogct = znum*(Act-1) + zct
    distsim(exogct,:,t) = 1.0

end do !zmin1ct
end do !zct
distsim(:,:,t) = distsim(:,:,t)/sum(distsim(:,:,t))

!do it for the firm-level stuff too
endogfirmpos(:,:) = 0
endogfirmpos(t,:) = endognum/2

!do this for the second period as well, because need two-period initialization here
t = 2
!distsim(:,:,t) = distsim(:,:,t-1)
Amin1ct = Asimpos(1)
Act = Asimpos(2)
do zmin1ct=1,1
  do zct=1,znum
    oldexogct = znum*(oldAct-1) + zct
    exogct = znum*(Act-1) + zct
    distsim(exogct,:,t) = distsim(oldexogct,:,t-1)

end do !zmin1ct
end do !zct

endogfirmpos(t,:) = endognum/2

!now do the endogenous portions for subsequent periods
do t=2,Tsim-1
  if (mod(t,250)==1) then
    write(*,*) "Doing period t = ",t
  end if

  Amin1ct=Asimpos(t-1)
  Act = Asimpos(t)
  Aprimect=Asimpos(t+1)


do zmin1ct=1,1
do zct=1,znum

  exogct = znum*(Act-1) + zct


  !loop over endog states
  do endogct=1,endognum
    if (distsim(exogct,endogct,t)>0.0) then


      !does the firm default?
      dfct = dfpolind(exogct,endogct)

      if (dfct==0) then

  		  kbpolct = kbpolind(exogct,endogct)

    	else if (dfct==1) then

    		kbpolct = 1

    	end if

      do zprimect=1,znum

        exogprimect = znum*(Aprimect-1) + zprimect

        distsim(exogprimect,kbpolct,t+1) = distsim(exogprimect,kbpolct,t+1) +&
          distsim(exogct,endogct,t)*pr_mat_z(zct,zprimect)

      end do !zprimect

  end if !distribution 0 check
  end do !endogct

end do !zct
end do !zmin1ct


distsim(:,:,t+1) = distsim(:,:,t+1)/sum(distsim(:,:,t+1))

!now, do the firm-level sim for the same period

do firmct=1,Nfirm

  !extract endog + exog stuff
  endogct = endogfirmpos(t,firmct)
  zct = zfirmpos(t,firmct)
  zmin1ct=zfirmpos(t-1,firmct)
  exogct = znum*(Act-1) + zct

  !does the firm default?
  dfct = dfpolind(exogct,endogct)

  !extract corresponding correct policies
  if (dfct==0) then

    kbpolct = kbpolind(exogct,endogct)

  else if (dfct==1) then

    kbpolct = 1

  end if

  kprimeval = endogval(kbpolct,1)
  bprimeval = endogval(kbpolct,2)
  kval = endogval(endogct,1)

  endogfirmpos(t+1,firmct) = kbpolct
  dfirmpos(t,firmct) = dfct
  kfirmpos(t,firmct) = endogind(endogct,1)
  bfirmpos(t,firmct) = endogind(endogct,2)
  pifirm(t,firmct) = pimat(exogct,endogct)+sigmapi*pifirmshocks(t,firmct) + sigmapi*gammapi*pifirmshocks(t-1,firmct)
  ifirm(t,firmct) = kprimeval - (1.0-delta)*kval
  kfirm(t,firmct) = endogval(endogct,1)
  bfirm(t,firmct) = endogval(endogct,2)
  !qfirm(t,firmct) = qmat(exogct,endogct)
  qfirm(t,firmct) = qmat(exogct,kbpolct)

  sprfirm(t,firmct) = 0.0
  if (qfirm(t,firmct)>0.0) then
    sprfirm(t,firmct) = (1.0/qfirm(t,firmct)) - (1.0+r)
  end if
  if (sprfirm(t,firmct)>0.5) then
    sprfirm(t,firmct) = 0.5
  end if

  pifcstfirm(t+1,firmct) = 0.0
  do Aprimect=1,Anum
    do zprimect=1,znum
      exogprimect = znum*(Aprimect-1) + zprimect
      pifcstfirm(t+1,firmct) = pifcstfirm(t+1,firmct) + &
        pimat(exogprimect,kbpolct)*pr_mat_exog_theta(exogct,exogprimect)
    end do
  end do

end do !firmct

end do !t

write(*,*) "Q, Spr mean vals = ",sum(qfirm)/dble(Tsim*Nfirm),sum(sprfirm)/dble(Tsim*Nfirm)

fefirm = pifirm -pifcstfirm

!now, fill in the aggregate variables
Ysim(:) = 0.0
Ksim(:) = 0.0
Isim(:) = 0.0
qsim(:) = 0.0
Bsim(:) = 0.0
dfsim(:) = 0.0
BtoKsim(:) = 0.0
irateheterosim(:,:) = 0.0
Dbheterosim(:,:) = 0.0
heterowgtvec(:) = 0.0
spreadsim(:) = 0.0
retbondsim(:) = 0.0
spreadtracksim(:,:) = 0.0
Nsim(:) = 0.0
sdMPKsim(:) = 0.0
dfYsim(:) = 0.0
dfBsim(:) = 0.0
dfIsim(:) = 0.0
issjunksim(:) =0.0
FEsim(:) = 0.0
Pisim(:) = 0.0
Pifsim(:) = 0.0

do t = 2,Tsim-1

  heterowgtvec(:) = 0.0

	Aval = Asim(t)

  Amin1ct=Asimpos(t-1)
  Act = Asimpos(t)


  issdenomval = 0.0
  fracissuing = 0.0
  EMPK = 0.0
  EMPK2 = 0.0
  junkdenomval = 0.0

  poskmass = 0.0
  poskprimemass = 0.0
  do zct=1,znum
    do zmin1ct=1,1

      exogct = znum*(Act-1) + zct
      zval = z0(zct)

      do endogct=1,endognum
        if (distsim(exogct,endogct,t)>0.0) then


        kval = endogval(endogct,1)
        bval = endogval(endogct,2)

        dfct = dfpolind(exogct,endogct)

        if (dfct==0) then

          kbpolct = kbpolind(exogct,endogct)

        else if (dfct==1) then

          kbpolct = 1

        end if

        kprimeval = endogval(kbpolct,1)
        bprimeval = endogval(kbpolct,2)
      	!qval = qmat(exogct,endogct)
        qval = qmat(exogct,kbpolct)
        nval = labor(Aval,zval,kval)
        yval = output(Aval,zval,kval)
        ival = kprimeval - (1.0-delta)*kval

        if (dfct==1) then
          yval = 0.0
        end if

        Ysim(t) = Ysim(t) + yval*distsim(exogct,endogct,t)
        Ksim(t) = Ksim(t) + kval*distsim(exogct,endogct,t)
      	Isim(t) = Isim(t) + ival*distsim(exogct,endogct,t)
      	qsim(t) = qsim(t) + qval*distsim(exogct,endogct,t)
      	Bsim(t) = Bsim(t) + bval*distsim(exogct,endogct,t)
      	dfsim(t) = dfsim(t) + dble(dfct)*distsim(exogct,endogct,t)
        if (kval>0.0) then
          BtoKsim(t) = BtoKsim(t) + (bval/kval)*distsim(exogct,endogct,t)
        end if
        spreadsim(t) = spreadsim(t) + (1.0/qval-(1.0+r))*distsim(exogct,endogct,t)
        !retbondsim(t) = retbondsim(t) + ((1.0-qsim(t-1))  + &
      !    dble(dfct)*(recovery(kval,bval,exogval(exogct,2),exogval(exogct,4))-1.0))&
      !      *distsim(exogct,endogct,t)
        retbondsim(t) = retbondsim(t) + ((1.0+ &
          dble(dfct)*(recovery(kval,bval,exogval(exogct,2),exogval(exogct,4))-1.0))/qsim(t-1))&
            *distsim(exogct,endogct,t)

        Nsim(t) = Nsim(t) + nval*distsim(exogct,endogct,t)
        EMPK = EMPK + log(yval/kval)*distsim(exogct,endogct,t)
        EMPK2 = EMPK2 + (log(yval/kval)**dble(2.0))*distsim(exogct,endogct,t)

        if (kval>0.0) then
        Pisim(t) = Pisim(t) + (pimat(exogct,endogct))*distsim(exogct,endogct,t)
        poskmass = poskmass + distsim(exogct,endogct,t)
        end if

        if (kprimeval>0.0) then
          do Aprimect=1,Anum
            do zprimect=1,znum
              exogprimect = znum*(Aprimect-1) + zprimect
              Pifsim(t+1) = Pifsim(t+1) + (pimat(exogprimect,kbpolct))*&
                pr_mat_exog_theta(exogct,exogprimect)*distsim(exogct,endogct,t)
            end do
          end do
          poskprimemass = poskprimemass + distsim(exogct,endogct,t)
        end if

            !hunt for correct location on spread grid and insert weight appropriately
            if (bprimeval > 0.0 ) then
              fracissuing = fracissuing + distsim(exogct,endogct,t)

            spreadval = (1.0/qval-(1.0+r))
            call hunt(spread0,spreadnum,spreadval,spreadind)

            if ((spreadind>=1).and.(spreadind<spreadnum)) then
              spreadwgt = (spreadval - spread0(spreadind))/&
                (spread0(spreadind+1)-spread0(spreadind))
            else if (spreadind<1) then
              spreadind = 1
              spreadwgt = 0.0
            else if (spreadind>=spreadnum) then
              spreadind = spreadnum-1
              spreadwgt = 1.0
            end if
            spreadtracksim(spreadind,t) = spreadtracksim(spreadind,t) + (1.0-spreadwgt)*distsim(exogct,endogct,t)
            spreadtracksim(spreadind+1,t) = spreadtracksim(spreadind+1,t) + spreadwgt*distsim(exogct,endogct,t)
            !b' > 0 selection
          end if



        !distributional selection
      end if
      end do !endogct

    end do !zmin1ct
  end do !zct

  spreadtracksim(:,t) = spreadtracksim(:,t) / fracissuing

  !Pisim(t) = Pisim(t)/poskmass
  !Pifsim(t) = Pifsim(t)/poskprimemass

  !process spread data in BAA spread
  sum40 = 0.0
  sum60 = 0.0
  found40 = 0
  found60 = 0

  do ct=2,spreadnum
    if (found40==0) then
      sum40 = sum40 + spreadtracksim(ct,t)/sum(spreadtracksim(2:spreadnum,t))
    end if

    if (found60==0) then
      sum60 = sum60 + spreadtracksim(ct,t)/sum(spreadtracksim(2:spreadnum,t))
    end if

    if ((sum40>0.4).and.(found40==0)) then
      found40 = ct
    end if

    if ((sum60>0.6).and.(found60==0)) then
      found60 = ct
    end if



  end do !ct

  spreadBAAsim(t) = sum(spread0(found40:found60)*spreadtracksim(found40:found60,t))/sum(spreadtracksim(found40:found60,t))



  sdMPKsim(t) = sqrt(EMPK2-EMPK**dble(2.0))

  spreadvec(:) = 0.0
  spr0ct =1
  do firmct=1,Nfirm
    if ((dfirmpos(t,firmct)==0).and.(bfirmpos(t+1,firmct)>0.0)) then

      spreadvec(firmct) = (1.0/qfirm(t,firmct))-(1.0+r)

    !end if !default conditional
    !if ((dfirmpos(t,firmct)==1).or.(spreadvec(firmct)<1.0e-4)) then
  else
      spr0ct=spr0ct + 1
    end if
  end do !firmct


  call heapsort(Nfirm,spreadvec)
  sprmedct = (spr0ct+Nfirm)/2




  !spreadsim(t) = spreadvec(sprmedct)

  !spreadsim(t) = spreadBAAsim(t)


end do !t
FEsim = Pisim - Pifsim

write(*,*) " "

if (profile_sigma_mode) then
  call profile_sigma_micro_moments(profile_m0,profile_mp,profile_mm)
  GMMmodelmom = profile_m0
else
!now, do the calculation of the GMM moments at the micro level

!first, compute the firm means
meani(:,:) = 0.0
Ni(:) = 0

do firmct=1,Nfirm
do t=Terg+1,Terg+Tfirm

  okobsflag = 1
  if (kfirm(t+2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t+1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-3,firmct)<=1.0e-4) then
    okobsflag = 0
  end if

  if (okobsflag==1) then

    Ni(firmct) = Ni(firmct) + 1


   !no firm dummies removed from the first 8 moments
   dumvec(1:8) = 0.0

    !fe/k_{t+1}
    dumvec(9) = (pifirm(t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)

    !pi/k_t
    dumvec(10) = pifirm(t,firmct)/kfirm(t,firmct)

    !i/k_t
    dumvec(11) = ifirm(t,firmct)/kfirm(t,firmct)

    !b'/k_t
    dumvec(12) = bfirm(t+1,firmct)/kfirm(t,firmct)

    !spread_t
    dumvec(13) = sprfirm(t,firmct)

    meani(firmct,:) = meani(firmct,:) + dumvec

!non-default selection
  end if

end do !t
end do !firmct

do firmct=1,Nfirm
meani(firmct,:) = meani(firmct,:)/dble(Ni(firmct))
end do !firmct

!then, compute the time means
meant(:,:) = 0.0
Nt(:) = 0

do firmct=1,Nfirm
do t=Terg+1,Terg+Tfirm

  okobsflag = 1
  if (kfirm(t+2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t+1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-3,firmct)<=1.0e-4) then
    okobsflag = 0
  end if

  if (okobsflag==1) then

    Ni(firmct) = Ni(firmct) + 1
    Nt(t-Terg) = Nt(t-Terg) + 1


    !fe/k_{t+1} - fe/k_t
    dumvec(1) = (pifirm(t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct) - &
      (pifirm(t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)

    !fe/k_{t+2} - fe/k_t
    dumvec(2) = (pifirm(t+2,firmct)-pifcstfirm(t+2,firmct))/kfirm(t+1,firmct) - &
      (pifirm(t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)

    !pi/k_t-1 - pi/k_t-2
    dumvec(3) = pifirm(t-1,firmct)/kfirm(t-1,firmct) - pifirm(t-2,firmct)/kfirm(t-2,firmct)

    !pi/k_t-1 - pi/k_t-3
    dumvec(4) = pifirm(t-1,firmct)/kfirm(t-1,firmct) - pifirm(t-3,firmct)/kfirm(t-3,firmct)

    !i/k_t-1 - i/k_t-2
    dumvec(5) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - ifirm(t-2,firmct)/kfirm(t-2,firmct)

    !i/k_t-1 - i/k_t-3
    dumvec(6) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - ifirm(t-3,firmct)/kfirm(t-3,firmct)

    !b'/k_t-1 - b'/k_t-2
    dumvec(7) = bfirm(t,firmct)/kfirm(t-1,firmct) - bfirm(t-1,firmct)/kfirm(t-2,firmct)

    !b'/k_t-1 - b'/k_t-3
    dumvec(8) = bfirm(t,firmct)/kfirm(t-1,firmct) - bfirm(t-2,firmct)/kfirm(t-3,firmct)

     !fe/k_{t+1}
     dumvec(9) = (pifirm(t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)

     !pi/k_t
     dumvec(10) = pifirm(t,firmct)/kfirm(t,firmct)

     !i/k_t
     dumvec(11) = ifirm(t,firmct)/kfirm(t,firmct)

     !b'/k_t
     dumvec(12) = bfirm(t+1,firmct)/kfirm(t,firmct)

     !spread_t
     dumvec(13) = sprfirm(t,firmct)

    dumvec = dumvec - meani(firmct,:)

    meant(t-Terg,:) = meant(t-Terg,:) + dumvec

!non-default selection
  end if

end do !t
end do !firmct

do t=1,Tfirm
meant(t,:) = meant(t,:)/dble(Nt(t))
end do !t

!then, compute the moments of the firm and time demeaned series
meanvec(:) = 0.0
covmat(:,:) = 0.0


ct = 0
do firmct=1,Nfirm
do t=Terg+1,Terg+Tfirm
  okobsflag = 1
  if (kfirm(t+2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t+1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-1,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-2,firmct)<=1.0e-4) then
    okobsflag = 0
  end if
  if (kfirm(t-3,firmct)<=1.0e-4) then
    okobsflag = 0
  end if

  if (okobsflag==1) then

    ct = ct + 1

    !fe/k_{t+1} - fe/k_t
    dumvec(1) = (pifirm(t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct) - &
      (pifirm(t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)

    !fe/k_{t+2} - fe/k_t
    dumvec(2) = (pifirm(t+2,firmct)-pifcstfirm(t+2,firmct))/kfirm(t+1,firmct) - &
      (pifirm(t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)

    !pi/k_t-1 - pi/k_t-2
    dumvec(3) = pifirm(t-1,firmct)/kfirm(t-1,firmct) - pifirm(t-2,firmct)/kfirm(t-2,firmct)

    !pi/k_t-1 - pi/k_t-3
    dumvec(4) = pifirm(t-1,firmct)/kfirm(t-1,firmct) - pifirm(t-3,firmct)/kfirm(t-3,firmct)

    !i/k_t-1 - i/k_t-2
    dumvec(5) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - ifirm(t-2,firmct)/kfirm(t-2,firmct)

    !i/k_t-1 - i/k_t-3
    dumvec(6) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - ifirm(t-3,firmct)/kfirm(t-3,firmct)

    !b'/k_t-1 - b'/k_t-2
    dumvec(7) = bfirm(t,firmct)/kfirm(t-1,firmct) - bfirm(t-1,firmct)/kfirm(t-2,firmct)

    !b'/k_t-1 - b'/k_t-3
    dumvec(8) = bfirm(t,firmct)/kfirm(t-1,firmct) - bfirm(t-2,firmct)/kfirm(t-3,firmct)

     !fe/k_{t+1}
     dumvec(9) = (pifirm(t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)

     !pi/k_t
     dumvec(10) = pifirm(t,firmct)/kfirm(t,firmct)

     !i/k_t
     dumvec(11) = ifirm(t,firmct)/kfirm(t,firmct)

     !b'/k_t
     dumvec(12) = bfirm(t+1,firmct)/kfirm(t,firmct)

     !spread_t
     dumvec(13) = sprfirm(t,firmct)



  dumvec = dumvec - meant(t-Terg,:) - meani(firmct,:)

  meanvec = meanvec + dumvec

  do ct1=1,numstats
    do ct2=1,numstats
      covmat(ct1,ct2) = covmat(ct1,ct2) + dumvec(ct1)*dumvec(ct2)
    end do !ct2
  end do !ct1

!non-default selection
  end if

end do !t
end do !firmct


!convert to moments rather than sums
meanvec = meanvec/dble(ct)
covmat = covmat/dble(ct)

!convert to centered moments
do ct1=1,numstats
  do ct2=1,numstats
    covmat(ct1,ct2) = covmat(ct1,ct2) -meanvec(ct1)*meanvec(ct2)
    end do !ct2
end do !ct1

!convert to correlations
do ct1=1,numstats
  do ct2=1,numstats
    corrmat(ct1,ct2) = covmat(ct1,ct2)/sqrt(covmat(ct1,ct1)*covmat(ct2,ct2))
    end do !ct2
end do !ct1


!if (1==0) then
write(*,*) "E vec = ",meanvec
write(*,*) " "
write(*,*) "Cov = "
write(*,*) covmat(1,:)
write(*,*) covmat(2,:)
write(*,*) covmat(3,:)
write(*,*) covmat(4,:)
write(*,*) covmat(5,:)
write(*,*) covmat(6,:)
write(*,*) covmat(7,:)
write(*,*) covmat(8,:)
write(*,*) covmat(9,:)
write(*,*) covmat(10,:)
write(*,*) covmat(11,:)
write(*,*) covmat(12,:)
write(*,*) covmat(13,:)


write(*,*) " "
write(*,*) "Corr = "
write(*,*) corrmat(1,:)
write(*,*) corrmat(2,:)
write(*,*) corrmat(3,:)
write(*,*) corrmat(4,:)
write(*,*) corrmat(5,:)
write(*,*) corrmat(6,:)
write(*,*) corrmat(7,:)
write(*,*) corrmat(8,:)
write(*,*) corrmat(9,:)
write(*,*) corrmat(10,:)
write(*,*) corrmat(11,:)
write(*,*) corrmat(12,:)
write(*,*) corrmat(13,:)


!end if


!now, insert the micro moments into the GMM vector
GMMmodelmom(:) = 0.0

!1 cov(d_su fe_{t+su}, d_sl pi_{t-1})
GMMmodelmom(1) = 0.25*(covmat(1,3) + covmat(1,4) + covmat(2,3) + covmat(2,4))

!2 cov(d_su fe_{t+su}, d_sl i_{t-1})
GMMmodelmom(2) = 0.25*(covmat(1,5) + covmat(1,6) + covmat(2,5) + covmat(2,6))

!3 cov(d_su fe_{t+su}, d_sl b'_{t-1})
GMMmodelmom(3) = 0.25*(covmat(1,7) + covmat(1,8) + covmat(2,7) + covmat(2,8))

!4 cov(fe_t+1 , fe_t+1)
GMMmodelmom(4) = covmat(9,9)

!5 cov(pi_t , pi_t)
GMMmodelmom(5) = covmat(10,10)

!6 cov(pi_t , i_t)
GMMmodelmom(6) = covmat(10,11)

!7 cov(pi_t , b'_t)
GMMmodelmom(7) = covmat(10,12)

!8 cov(pi_t , spr_t)
GMMmodelmom(8) = covmat(10,13)

!9 cov(i_t , i_t)
GMMmodelmom(9) = covmat(11,11)

!10 cov(i_t , b'_t)
GMMmodelmom(10) = covmat(11,12)

!11 cov(i_t , spr_t)
GMMmodelmom(11) = covmat(11,13)

!12 cov(b'_t , b'_t)
GMMmodelmom(12) = covmat(12,12)

!13 cov(b'_t , spr_t)
GMMmodelmom(13) = covmat(12,13)

!14 cov(spr_t , spr_t)
GMMmodelmom(14) = covmat(13,13)

end if

!now, do the calculations of the GMM moments at the macro level
!dumvecmacro(numstatsmacro),meanvecmacro(numstatsmacro),&
!covmatmacro(numstatsmacro,numstatsmacro),corrmatmacro(numstatsmacro,numstatsmacro)

meanvecmacro(:) = 0.0
covmatmacro(:,:) = 0.0
do t=Terg,Tsim-2
    dumvecmacro = (/ spreadsim(t), &
                  dfsim(t),&
                  log(Ysim(t))-log(Ysim(t-1))/)

    meanvecmacro = meanvecmacro + dumvecmacro
    do ct1=1,numstatsmacro
      do ct2=1,numstatsmacro
        covmatmacro(ct1,ct2) = covmatmacro(ct1,ct2)  + dumvecmacro(ct1)*dumvecmacro(ct2)
      end do !ct2
    end do !ct1

end do !t

meanvecmacro = meanvecmacro / dble(Tsim-2 - Terg + 1)
covmatmacro = covmatmacro / dble(Tsim-2 - Terg + 1)

!convert to centered moments
do ct1=1,numstatsmacro
  do ct2=1,numstatsmacro
    covmatmacro(ct1,ct2) = covmatmacro(ct1,ct2) -meanvecmacro(ct1)*meanvecmacro(ct2)
    end do !ct2
end do !ct1

!convert to correlations
do ct1=1,numstatsmacro
  do ct2=1,numstatsmacro
    corrmatmacro(ct1,ct2) = covmatmacro(ct1,ct2)/sqrt(covmatmacro(ct1,ct1)*covmatmacro(ct2,ct2))
    end do !ct2
end do !ct1

!if (1==0) then
write(*,*) " "
write(*,*) "E vec (macro) = ",meanvecmacro
write(*,*) " "

write(*,*) "Cov (macro) = "
write(*,*) covmatmacro(1,:)
write(*,*) covmatmacro(2,:)
write(*,*) covmatmacro(3,:)
write(*,*) " "

write(*,*) "Corr (macro) = "
write(*,*) corrmatmacro(1,:)
write(*,*) corrmatmacro(2,:)
write(*,*) corrmatmacro(3,:)
write(*,*) " "

!end if


!insert moments into the GMM vector at the macro level
GMMmodelmom(15) = meanvecmacro(1) !spread
GMMmodelmom(16) = meanvecmacro(2) !default
GMMmodelmom(17) = covmatmacro(3,3) !cov(D GDP,D GDP)

if (profile_sigma_mode) then
  profile_m0(15:17) = GMMmodelmom(15:17)
  profile_mp(15:17) = GMMmodelmom(15:17)
  profile_mm(15:17) = GMMmodelmom(15:17)
  GMMmodelmom = profile_m0
end if

!this set of code computes the objective according to the optimal weighting matrix
do ct=1,GMMmom_small
  momdiffvec(ct) = GMMdatamom(GMMmomselect(ct))-GMMmodelmom(GMMmomselect(ct))
end do !ct

do ct1=1,GMMmom_small
  do ct2=1,GMMmom_small
    weightmat(ct1,ct2) = GMMdatawgt(GMMmomselect(ct1),GMMmomselect(ct2))
  end do !ct2
end do !ct1

!momdiffvec = GMMdatamom - GMMmodelmom
quadformvec = matmul(weightmat,momdiffvec)
GMMobjval = dot_product(momdiffvec,quadformvec)

write(*,*) "GMM calcs: Data, Model"
do ct=1,GMMmom
write(*,*) GMMdatamom(ct),GMMmodelmom(ct)
end do !ct
write(*,*) " "
write(*,*) "GMM obj = ",GMMobjval
write(*,*) " "

end subroutine simulate

subroutine profile_sigma_micro_moments(profile_m0,profile_mp,profile_mm)
use params
use base_lib
implicit none

double precision, intent(out) :: profile_m0(GMMmom),profile_mp(GMMmom),profile_mm(GMMmom)
integer :: ct,t,firmct,lane,ct1,ct2,okobsflag,Nt(Tfirm),Ni(Nfirm)
double precision, allocatable :: pifirm_profile(:,:,:)
double precision :: dumvec_profile(3,numstats),meanvec_profile(3,numstats),&
  covmat_profile(3,numstats,numstats),meani_profile(3,Nfirm,numstats),&
  meant_profile(3,Tfirm,numstats),moments_profile(3,GMMmom)

! Build the three measured-profit lanes with the released operation order.
! The policy, state, shock, forecast-profit, and aggregate paths are shared.
allocate(pifirm_profile(3,Tsim,Nfirm))
pifirm_profile(:,:,:) = 0.0d0
do firmct=1,Nfirm
  do t=2,Tsim-1
    pifirm_profile(1,t,firmct) = pifirm(t,firmct)
    pifirm_profile(2,t,firmct) = pifirm(t,firmct) + &
      pifirmshocks(t,firmct) + gammapi*pifirmshocks(t-1,firmct)
    pifirm_profile(3,t,firmct) = pifirm(t,firmct) - &
      pifirmshocks(t,firmct) - gammapi*pifirmshocks(t-1,firmct)
  end do
end do

meani_profile(:,:,:) = 0.0d0
Ni(:) = 0
do firmct=1,Nfirm
  do t=Terg+1,Terg+Tfirm
    okobsflag = 1
    if (kfirm(t+2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t+1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-3,firmct)<=1.0e-4) okobsflag = 0
    if (okobsflag==1) then
      Ni(firmct) = Ni(firmct) + 1
      do lane=1,3
        dumvec_profile(lane,1:8) = 0.0d0
        dumvec_profile(lane,9) = &
          (pifirm_profile(lane,t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)
        dumvec_profile(lane,10) = pifirm_profile(lane,t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,11) = ifirm(t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,12) = bfirm(t+1,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,13) = sprfirm(t,firmct)
        meani_profile(lane,firmct,:) = meani_profile(lane,firmct,:) + &
          dumvec_profile(lane,:)
      end do
    end if
  end do
end do
do lane=1,3
  do firmct=1,Nfirm
    meani_profile(lane,firmct,:) = meani_profile(lane,firmct,:)/dble(Ni(firmct))
  end do
end do

meant_profile(:,:,:) = 0.0d0
Nt(:) = 0
do firmct=1,Nfirm
  do t=Terg+1,Terg+Tfirm
    okobsflag = 1
    if (kfirm(t+2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t+1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-3,firmct)<=1.0e-4) okobsflag = 0
    if (okobsflag==1) then
      Nt(t-Terg) = Nt(t-Terg) + 1
      do lane=1,3
        dumvec_profile(lane,1) = &
          (pifirm_profile(lane,t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct) - &
          (pifirm_profile(lane,t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)
        dumvec_profile(lane,2) = &
          (pifirm_profile(lane,t+2,firmct)-pifcstfirm(t+2,firmct))/kfirm(t+1,firmct) - &
          (pifirm_profile(lane,t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)
        dumvec_profile(lane,3) = pifirm_profile(lane,t-1,firmct)/kfirm(t-1,firmct) - &
          pifirm_profile(lane,t-2,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,4) = pifirm_profile(lane,t-1,firmct)/kfirm(t-1,firmct) - &
          pifirm_profile(lane,t-3,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,5) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - &
          ifirm(t-2,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,6) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - &
          ifirm(t-3,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,7) = bfirm(t,firmct)/kfirm(t-1,firmct) - &
          bfirm(t-1,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,8) = bfirm(t,firmct)/kfirm(t-1,firmct) - &
          bfirm(t-2,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,9) = &
          (pifirm_profile(lane,t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)
        dumvec_profile(lane,10) = pifirm_profile(lane,t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,11) = ifirm(t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,12) = bfirm(t+1,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,13) = sprfirm(t,firmct)
        dumvec_profile(lane,:) = dumvec_profile(lane,:) - &
          meani_profile(lane,firmct,:)
        meant_profile(lane,t-Terg,:) = meant_profile(lane,t-Terg,:) + &
          dumvec_profile(lane,:)
      end do
    end if
  end do
end do
do lane=1,3
  do t=1,Tfirm
    meant_profile(lane,t,:) = meant_profile(lane,t,:)/dble(Nt(t))
  end do
end do

meanvec_profile(:,:) = 0.0d0
covmat_profile(:,:,:) = 0.0d0
ct = 0
do firmct=1,Nfirm
  do t=Terg+1,Terg+Tfirm
    okobsflag = 1
    if (kfirm(t+2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t+1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-1,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-2,firmct)<=1.0e-4) okobsflag = 0
    if (kfirm(t-3,firmct)<=1.0e-4) okobsflag = 0
    if (okobsflag==1) then
      ct = ct + 1
      do lane=1,3
        dumvec_profile(lane,1) = &
          (pifirm_profile(lane,t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct) - &
          (pifirm_profile(lane,t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)
        dumvec_profile(lane,2) = &
          (pifirm_profile(lane,t+2,firmct)-pifcstfirm(t+2,firmct))/kfirm(t+1,firmct) - &
          (pifirm_profile(lane,t,firmct)-pifcstfirm(t,firmct))/kfirm(t-1,firmct)
        dumvec_profile(lane,3) = pifirm_profile(lane,t-1,firmct)/kfirm(t-1,firmct) - &
          pifirm_profile(lane,t-2,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,4) = pifirm_profile(lane,t-1,firmct)/kfirm(t-1,firmct) - &
          pifirm_profile(lane,t-3,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,5) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - &
          ifirm(t-2,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,6) = ifirm(t-1,firmct)/kfirm(t-1,firmct) - &
          ifirm(t-3,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,7) = bfirm(t,firmct)/kfirm(t-1,firmct) - &
          bfirm(t-1,firmct)/kfirm(t-2,firmct)
        dumvec_profile(lane,8) = bfirm(t,firmct)/kfirm(t-1,firmct) - &
          bfirm(t-2,firmct)/kfirm(t-3,firmct)
        dumvec_profile(lane,9) = &
          (pifirm_profile(lane,t+1,firmct)-pifcstfirm(t+1,firmct))/kfirm(t,firmct)
        dumvec_profile(lane,10) = pifirm_profile(lane,t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,11) = ifirm(t,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,12) = bfirm(t+1,firmct)/kfirm(t,firmct)
        dumvec_profile(lane,13) = sprfirm(t,firmct)
        dumvec_profile(lane,:) = dumvec_profile(lane,:) - &
          meant_profile(lane,t-Terg,:) - meani_profile(lane,firmct,:)
        meanvec_profile(lane,:) = meanvec_profile(lane,:) + dumvec_profile(lane,:)
        do ct1=1,numstats
          do ct2=1,numstats
            covmat_profile(lane,ct1,ct2) = covmat_profile(lane,ct1,ct2) + &
              dumvec_profile(lane,ct1)*dumvec_profile(lane,ct2)
          end do
        end do
      end do
    end if
  end do
end do

do lane=1,3
  meanvec_profile(lane,:) = meanvec_profile(lane,:)/dble(ct)
  covmat_profile(lane,:,:) = covmat_profile(lane,:,:)/dble(ct)
  do ct1=1,numstats
    do ct2=1,numstats
      covmat_profile(lane,ct1,ct2) = covmat_profile(lane,ct1,ct2) - &
        meanvec_profile(lane,ct1)*meanvec_profile(lane,ct2)
    end do
  end do
  moments_profile(lane,:) = 0.0d0
  moments_profile(lane,1) = 0.25d0*(covmat_profile(lane,1,3) + &
    covmat_profile(lane,1,4) + covmat_profile(lane,2,3) + covmat_profile(lane,2,4))
  moments_profile(lane,2) = 0.25d0*(covmat_profile(lane,1,5) + &
    covmat_profile(lane,1,6) + covmat_profile(lane,2,5) + covmat_profile(lane,2,6))
  moments_profile(lane,3) = 0.25d0*(covmat_profile(lane,1,7) + &
    covmat_profile(lane,1,8) + covmat_profile(lane,2,7) + covmat_profile(lane,2,8))
  moments_profile(lane,4) = covmat_profile(lane,9,9)
  moments_profile(lane,5) = covmat_profile(lane,10,10)
  moments_profile(lane,6) = covmat_profile(lane,10,11)
  moments_profile(lane,7) = covmat_profile(lane,10,12)
  moments_profile(lane,8) = covmat_profile(lane,10,13)
  moments_profile(lane,9) = covmat_profile(lane,11,11)
  moments_profile(lane,10) = covmat_profile(lane,11,12)
  moments_profile(lane,11) = covmat_profile(lane,11,13)
  moments_profile(lane,12) = covmat_profile(lane,12,12)
  moments_profile(lane,13) = covmat_profile(lane,12,13)
  moments_profile(lane,14) = covmat_profile(lane,13,13)
end do

profile_m0 = moments_profile(1,:)
profile_mp = moments_profile(2,:)
profile_mm = moments_profile(3,:)
deallocate(pifirm_profile)
end subroutine profile_sigma_micro_moments

subroutine distribution()
use params
use base_lib
use omp_lib
implicit none

integer :: exogct,endogct,zct,Act,zprimect,Aprimect,exogprimect,kbpolct,&
	kprimect,bprimect,distct,ind,kct,bct
double precision :: bval,kval,kprimeval,disterr,dumkvec(knum),dumbvec(bnum),wgt

write(*,*) "##### Solving for ergodic distribution "


distold(:,:) = 1.0
distold = distold/sum(distold)


do distct=1,distmaxit

!loop over states
!NOTE: tried parallelization here, no gain because of overhead
do exogct=1,exognum
do endogct=1,endognum
if (distold(exogct,endogct)>0.0) then

	!what are the policies & states?
	kval = endogval(endogct,1)
	kbpolct=kbpolind(exogct,endogct)
	Act = exogind(exogct,2)
	zct = exogind(exogct,4)

	!if you don't default
	if (dfpolind(exogct,endogct)==0) then

		do Aprimect=1,Anum
			do zprimect=1,znum
				exogprimect = znum*(Aprimect-1) + zprimect
				dist(exogprimect,kbpolct) = dist(exogprimect,kbpolct) + &
					pr_mat_exog(exogct,exogprimect)*distold(exogct,endogct)
			end do !zprimect
		end do !Aprimect

	end if

	!if you do default
	if (dfpolind(exogct,endogct)==1) then
		kprimect = 1
		bprimect = 1

		kbpolct = bnum*(kprimect-1) + bprimect

		do Aprimect=1,Anum
			do zprimect=1,znum
				exogprimect = znum*(Aprimect-1) + zprimect
				dist(exogprimect,kbpolct) = dist(exogprimect,kbpolct) + &
					pr_mat_exog(exogct,exogprimect)*distold(exogct,endogct)

			end do !zprimect
		end do !Aprimect

	end if


end if
end do !endogct
end do !exogct


disterr = maxval(dist-distold)

if (mod(distct,10)==1) then
write(*,*) "distct = ",distct," err = ",disterr
end if

if (disterr<disttol) then
	distFLAG = 0
	exit
end if

dist = dist/sum(dist)
distold = dist
dist(:,:) = 0.0

end do !distct

distk(:) = 0.0
distb(:) = 0.0
do kct=1,knum
	do bct=1,bnum
		endogct=bnum*(kct-1) + bct
		distk(kct) = distk(kct) + sum(dist(:,endogct))
		distb(bct) = distb(bct) + sum(dist(:,endogct))
	end do !bct
end do !kct

write(*,*) " "

bunchFLAG =0
if ((distb(bnum)+distk(2)+distk(knum))>0.05) then
  bunchFLAG=1
end if


end subroutine distribution

subroutine expected_value_table(Vin,EVout)
  use params
  use omp_lib
  implicit none

  double precision, intent(in) :: Vin(exognum,endognum)
  double precision, intent(out) :: EVout(exognum,endognum)
  integer :: exogct,kbpolct,Act,zct,Aprimect,zprimect,exogprimect

  !$omp parallel do collapse(2) default(shared) &
  !$omp& private(exogct,kbpolct,Act,zct,Aprimect,zprimect,exogprimect)
  do kbpolct=1,endognum
    do exogct=1,exognum
      Act = exogind(exogct,2)
      zct = exogind(exogct,4)
      EVout(exogct,kbpolct) = 0.0
      do Aprimect=1,Anum
        do zprimect=1,znum
          exogprimect = znum*(Aprimect-1) + zprimect
          if (doRATfirms==0) then
            EVout(exogct,kbpolct) = EVout(exogct,kbpolct) + beta * &
              pr_mat_exog_theta(exogct,exogprimect) * &
              Vin(exogprimect,kbpolct)
          else if (doRATfirms==1) then
            EVout(exogct,kbpolct) = EVout(exogct,kbpolct) + beta * &
              pr_mat_exog(exogct,exogprimect) * Vin(exogprimect,kbpolct)
          end if
        end do
      end do
    end do
  end do
  !$omp end parallel do
end subroutine expected_value_table

subroutine default_value_reference_order(Vin,Vdout)
  use params
  use omp_lib
  implicit none

  double precision, intent(in) :: Vin(exognum,endognum)
  double precision, intent(out) :: Vdout(exognum)
  integer :: exogct,exogprimect,Act,zct,Aprimect,zprimect

  do exogct=1,exognum
    Act = exogind(exogct,2)
    zct = exogind(exogct,4)
    Vdout(exogct) = 0.0
    if (doRATfirms==0) then
      do Aprimect=1,Anum
        do zprimect=1,znum
          exogprimect = znum*(Aprimect-1) + zprimect
          Vdout(exogct) = Vdout(exogct) + &
            beta*pr_mat_exog_theta(exogct,exogprimect)*Vin(exogprimect,1)
        end do
      end do
    else if (doRATfirms==1) then
      do Aprimect=1,Anum
        do zprimect=1,znum
          exogprimect = znum*(Aprimect-1) + zprimect
          Vdout(exogct) = Vdout(exogct) + &
            beta*pr_mat_exog(exogct,exogprimect)*Vin(exogprimect,1)
        end do
      end do
    end if
  end do
end subroutine default_value_reference_order

subroutine politeration()
	use params
	use base_lib
	use omp_lib
	use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
	implicit none

	integer :: zmin1ct,zct,Amin1ct,Act,polct,vfct,kbpolct,zprimect,Aprimect,&
		exogprimect,endogct,exogct,dfct,accelct,bprimect,kprimect,bct,kct,&
		topct,ntop,topind(endognum),active_count,active_pos,push_count,&
		active_index(endognum),push_index(endognum),active_updates,active_min,active_max
	integer(kind=8) :: active_total
	double precision :: zmin1val,zval,Amin1val,Aval,dfprimeval,kprimeval,bprimeval,&
		kval,bval,RHSvec(endognum),RHSndvec(endognum),qval,dval,polerr,vferr,dferr,&
	    ACval,ival,yval,dumbvec(bnum),dumkvec(knum),Aprimeval,zprimeval,dfmeanerr,&
	    minshift,maxshift,RHSdivvec(endognum),approxmax
	double precision, allocatable :: EVtheta(:,:),EVthetaT(:,:),qmatT(:,:),&
		currentdiv(:,:),ynetcache(:,:,:),accache(:,:),currentbase(:,:,:),&
		actionkbase(:,:),debtcost(:,:),actioncost(:,:),actionb(:)
	double precision, allocatable :: lazy_prev_ev(:,:),lazy_gap(:,:),lazy_range(:)
	logical, allocatable :: actioninvalid(:)
	logical :: active_action(endognum)
	double precision, parameter :: tie_guard = 1.0d-12
	integer :: splitct,splitlo,splithi,splitmid,candidatect,actionsplit(endognum),&
		posbincount(0:bnum),negbincount(0:bnum),&
		poshead(0:bnum),postail(0:bnum),neghead(0:bnum),negtail(0:bnum),&
		posnext(endognum),negnext(endognum),bestaction
	integer(kind=8) :: envelope_candidates
	double precision :: dividendslope,actionintercept,rawdividend,bestapprox,&
		regimeapprox,currentdebtshift,&
		posconstant,negconstant,posbinval(0:bnum),negbinval(0:bnum),&
		posactionval(endognum),negactionval(endognum),&
		sparsevalue(endognum),sparsedividend(endognum),bestvalue,bestdividend
	logical :: sparse_nonfinite
	double precision :: lazy_secondbest,lazy_margin
	double precision, parameter :: lazy_roundoff_guard = 1.0d-14
	logical :: lazy_valid,lazy_block
	integer(kind=8) :: lazy_blocks,lazy_total_blocks,lazy_certified_states,lazy_total_states

	allocate(EVtheta(exognum,endognum),EVthetaT(endognum,exognum),&
		qmatT(endognum,exognum),currentdiv(exognum,endognum),&
		ynetcache(Anum,znum,knum),accache(knum,knum),&
		currentbase(Anum,znum,endognum),actionkbase(knum,knum),&
		debtcost(bnum,bnum),actioncost(endognum,endognum),&
		actionb(endognum),actioninvalid(endognum))
	allocate(lazy_prev_ev(endognum,exognum),lazy_gap(exognum,endognum),&
		lazy_range(exognum))

	! These terms are invariant within one model evaluation.  Retain the
	! released functions for their exact floating-point values, then reuse the
	! stored values in the approximate all-action scan.  The winning tie set is
	! still recomputed below with the released dividend function.
	do Act=1,Anum
		do zct=1,znum
			do kct=1,knum
				ynetcache(Act,zct,kct) = outputnetlabor(A0(Act),z0(zct),k0(kct))
			end do
		end do
	end do
	do kct=1,knum
		do kprimect=1,knum
			accache(kct,kprimect) = AC(k0(kct),k0(kprimect))
			actionkbase(kct,kprimect) = -(1.0-tau)*accache(kct,kprimect) &
				- invprice*k0(kprimect)
		end do
	end do
	do bct=1,bnum
		do bprimect=1,bnum
			debtcost(bct,bprimect) = 0.0
			if (b0(bprimect)>b0(bct)) then
				debtcost(bct,bprimect) = -(etab/2.0)*(b0(bprimect)-b0(bct))**dble(2.0)
			end if
		end do
	end do
	do kbpolct=1,endognum
		actionb(kbpolct) = endogval(kbpolct,2)
		actioninvalid(kbpolct) = endogind(kbpolct,1)==1
		do endogct=1,endognum
			kct = endogind(endogct,1)
			bct = endogind(endogct,2)
			actioncost(kbpolct,endogct) = &
				actionkbase(kct,endogind(kbpolct,1)) &
				+ debtcost(bct,endogind(kbpolct,2))
		end do
	end do
	do Act=1,Anum
		do zct=1,znum
			do endogct=1,endognum
				kct = endogind(endogct,1)
				bct = endogind(endogct,2)
				kval = endogval(endogct,1)
				bval = endogval(endogct,2)
				currentbase(Act,zct,endogct) = &
					(1.0-tau)*(ynetcache(Act,zct,kct)-phi) &
					+ invprice*(1.0-delta)*kval &
					+ tau*(r*bval+delta*kval)-bval
			end do
		end do
	end do

active_total = 0_8
active_updates = 0
active_min = endognum
active_max = 0
lazy_valid = .false.
lazy_blocks = 0_8
lazy_total_blocks = 0_8
lazy_certified_states = 0_8
lazy_total_states = 0_8
lazy_prev_ev(:,:) = 0.0d0
lazy_gap(:,:) = 0.0d0
lazy_range(:) = huge(1.0d0)

!do you actually do the VFI, or just read in policies + values?
if (doVFI==1) then

	write(*,*) "##### Doing policy iteration "
	write(*,*) " "

	dfpolindold(:,:) = 0
  kbpolindold(:,:) = endognum/2
	Vthetaold(:,:) = 0.0

	!outer loop - on the default states
	do dfct=1,dfmaxit
		lazy_valid = .false.

		write(*,*) "Default outer loop # ",dfct

		!given dfpolindold, compute the implied bond prices qmat
		qmat(:,:) = 1.0
		!$omp parallel private(exogct,kbpolct,Act,zct,kprimeval,bprimeval,&
		!$omp& Aprimect,zprimect,exogprimect,dfprimeval,Aprimeval,zprimeval)
		!$omp do collapse(2)
		do exogct=1,exognum
			do kbpolct=1,endognum

				!what are current exog states?
				Act = exogind(exogct,2)
				zct = exogind(exogct,4)

				!what is the optimal policy?
				kprimeval = endogval(kbpolct,1)
				bprimeval = endogval(kbpolct,2)

				if (bprimeval>0.0) then
					!loop over exog outcomes tomorrow
					do Aprimect=1,Anum
						do zprimect=1,znum
							exogprimect = znum*(Aprimect-1) + zprimect

              Aprimeval = exogval(exogprimect,2)
              zprimeval = exogval(exogprimect,4)

							dfprimeval = dble(dfpolindold(exogprimect,kbpolct))

              if (doRATprice==0) then
							qmat(exogct,kbpolct) = qmat(exogct,kbpolct) + &
								pr_mat_exog_theta(exogct,exogprimect)*&
								dfprimeval * (recovery(kprimeval,bprimeval,Aprimeval,zprimeval)-1.0)
              else if (doRATprice==1) then
                qmat(exogct,kbpolct) = qmat(exogct,kbpolct) + &
                  pr_mat_exog(exogct,exogprimect)*&
                  dfprimeval * (recovery(kprimeval,bprimeval,Aprimeval,zprimeval)-1.0)
              endif

						end do !zprimect
					end do !Aprimect

					qmat(exogct,kbpolct) = qmat(exogct,kbpolct)/(1.0+r)


				else if (bprimeval<=0.0) then

					qmat(exogct,kbpolct) = 1.0
					qmat(exogct,kbpolct) = qmat(exogct,kbpolct)/(1.0+r)

				end if


			end do !kbpolct
		end do !exogct
		!$omp end do nowait
		!$omp end parallel
		qmatT = transpose(qmat)

		write(*,*) " "

		!now, given the default states & bond pricing, determine new policies plus values
		do polct=1,polmaxit

			! The current policy and qmat do not change during the 100 value
			! push-forwards.  The released code recomputed this same dividend in
			! every pass; compute it once without changing its operation order.
			!$omp parallel do collapse(2) default(shared) &
			!$omp& private(exogct,endogct,Act,zct,kval,bval,kbpolct,kprimeval,bprimeval,qval,dval)
			do exogct=1,exognum
				do endogct=1,endognum
					Act = exogind(exogct,2)
					zct = exogind(exogct,4)
					kval = endogval(endogct,1)
					bval = endogval(endogct,2)
					kbpolct = kbpolindold(exogct,endogct)
					kprimeval = endogval(kbpolct,1)
					bprimeval = endogval(kbpolct,2)
					qval = qmat(exogct,kbpolct)
					dval = dividend(A0(Act),z0(zct),kval,kprimeval,qval,bprimeval,bval)
					currentdiv(exogct,endogct) = dval
					if (dval>=0.0) then
						divposmat(exogct,endogct) = 1.0
					else
						divposmat(exogct,endogct) = 0.0
					end if
				end do
			end do
			!$omp end parallel do

      !if (1==0) then
			! The fixed policy can read only endogenous states selected by a
			! continuing state. Default reads reset state one. Restrict the first
			! 99 Howard pushes to that live union, then reconstruct every state in
			! the final push before the complete policy improvement.
			active_action(:) = .false.
			active_action(1) = .true.
			do exogct=1,exognum
				do endogct=1,endognum
					if (dfpolindold(exogct,endogct)==0) then
						active_action(kbpolindold(exogct,endogct)) = .true.
					end if
				end do
			end do
			active_count = 0
			do endogct=1,endognum
				if (active_action(endogct)) then
					active_count = active_count+1
					active_index(active_count) = endogct
				end if
			end do
			active_updates = active_updates+1
			active_total = active_total+int(active_count,kind=8)
			active_min = min(active_min,active_count)
			active_max = max(active_max,active_count)

			!push the value function forward accelmaxit times
			do accelct=1,accelmaxit

				call default_value_reference_order(Vthetaold,Vdtheta)
				if (accelct<accelmaxit) then
					push_count = active_count
					push_index(1:push_count) = active_index(1:push_count)
				else
					push_count = endognum
					do endogct=1,endognum
						push_index(endogct) = endogct
					end do
				end if

				!$omp parallel private(active_pos,exogct,endogct,kbpolct,dval,&
				!$omp& Aprimect,zprimect,exogprimect)
				!$omp do collapse(2)
				do active_pos=1,push_count
					do exogct=1,exognum
						endogct = push_index(active_pos)
						kbpolct = kbpolindold(exogct,endogct)
						dval = currentdiv(exogct,endogct)

						if (dfpolindold(exogct,endogct)==0) then
							Vndtheta(exogct,endogct) = dval
							do Aprimect=1,Anum
								do zprimect=1,znum
									exogprimect = znum*(Aprimect-1)+zprimect
									if (doRATfirms==0) then
										Vndtheta(exogct,endogct) = Vndtheta(exogct,endogct) + beta * &
											pr_mat_exog_theta(exogct,exogprimect) * &
											Vthetaold(exogprimect,kbpolct)
									else if (doRATfirms==1) then
										Vndtheta(exogct,endogct) = Vndtheta(exogct,endogct) + beta * &
											pr_mat_exog(exogct,exogprimect) * &
											Vthetaold(exogprimect,kbpolct)
									end if
								end do
							end do
							Vtheta(exogct,endogct) = Vndtheta(exogct,endogct)
						else
							Vtheta(exogct,endogct) = Vdtheta(exogct)
						end if
					end do
				end do
				!$omp end do nowait
				!$omp end parallel

				Vthetaold = Vtheta

			end do !accelct
	    !end if

			call expected_value_table(Vthetaold,EVtheta)
			EVthetaT = transpose(EVtheta)
			call default_value_reference_order(Vthetaold,Vdtheta)
			if (lazy_valid) then
				do exogct=1,exognum
					lazy_range(exogct) = &
						maxval(EVthetaT(bnum+1:endognum,exogct)-&
							lazy_prev_ev(bnum+1:endognum,exogct))-&
						minval(EVthetaT(bnum+1:endognum,exogct)-&
							lazy_prev_ev(bnum+1:endognum,exogct))
				end do
			else
				lazy_range(:) = huge(1.0d0)
			end if

			! For etab=0, inherited debt enters every action's raw dividend with
			! the same slope tau*r-1.  Conditional on whether equity is issued,
			! action rankings are therefore independent of current debt.  Bin each
			! action by the last debt node with nonnegative raw dividends, retain
			! each bin's numerical tie set in both payout regimes, and recompute
			! every retained candidate with the released dividend and successor-add
			! order.  This is exact in real arithmetic; finite-precision parity is
			! calibrated separately against the operational authority.
			envelope_candidates = 0
			dividendslope = tau*r-1.0d0
			if (dividendslope>=0.0d0) &
				error stop 'debt-envelope monotonicity requires tau*r-1 < 0'
			!$omp parallel do collapse(2) default(shared) reduction(+:envelope_candidates) &
			!$omp& private(exogct,kct,Act,zct,kval,kbpolct,kprimect,bprimeval,qval,&
			!$omp& actionintercept,rawdividend,splitct,splitlo,splithi,splitmid,&
			!$omp& posconstant,negconstant,bestapprox,regimeapprox,currentdebtshift,&
			!$omp& posbinval,negbinval,actionsplit,posactionval,negactionval,&
			!$omp& posbincount,negbincount,poshead,postail,neghead,negtail,&
			!$omp& posnext,negnext,bct,endogct,candidatect,topind,ntop,topct,dval,&
			!$omp& sparsevalue,sparsedividend,bestvalue,bestdividend,bestaction,&
			!$omp& sparse_nonfinite,RHSndvec,RHSdivvec,Aprimect,zprimect,exogprimect,bval,&
			!$omp& lazy_block,lazy_secondbest,lazy_margin) &
			!$omp& reduction(+:lazy_blocks,lazy_total_blocks,lazy_certified_states,lazy_total_states)
			do exogct=1,exognum
				do kct=1,knum
					Act = exogind(exogct,2)
					zct = exogind(exogct,4)
					kval = k0(kct)
					lazy_total_blocks = lazy_total_blocks+1_8
					lazy_total_states = lazy_total_states+int(bnum,kind=8)
					if (lazy_valid) lazy_certified_states = lazy_certified_states+int(count(&
						lazy_gap(exogct,bnum*(kct-1)+1:bnum*kct)>&
						lazy_range(exogct)+lazy_roundoff_guard),kind=8)
					lazy_block = lazy_valid .and. all(&
						lazy_gap(exogct,bnum*(kct-1)+1:bnum*kct)>&
						lazy_range(exogct)+lazy_roundoff_guard)
					if (lazy_block) then
						lazy_blocks = lazy_blocks+1_8
						do bct=1,bnum
							endogct = bnum*(kct-1)+bct
							bval = b0(bct)
							kbpolct = kbpolindold(exogct,endogct)
							kprimect = endogind(kbpolct,1)
							bprimeval = endogval(kbpolct,2)
							qval = qmat(exogct,kbpolct)
							dval = dividend(A0(Act),z0(zct),kval,k0(kprimect),&
								qval,bprimeval,bval)
							bestvalue = dval
							do Aprimect=1,Anum
								do zprimect=1,znum
									exogprimect = znum*(Aprimect-1)+zprimect
									if (doRATfirms==0) then
										bestvalue = bestvalue+beta*&
											pr_mat_exog_theta(exogct,exogprimect)*&
											Vthetaold(exogprimect,kbpolct)
									else if (doRATfirms==1) then
										bestvalue = bestvalue+beta*&
											pr_mat_exog(exogct,exogprimect)*&
											Vthetaold(exogprimect,kbpolct)
									end if
								end do
							end do
							kbpolind(exogct,endogct) = kbpolct
							Vndtheta(exogct,endogct) = bestvalue
							if (dval>=0.0d0) then
								divposmat(exogct,endogct) = 1.0d0
							else
								divposmat(exogct,endogct) = 0.0d0
							end if
							if (bestvalue>=Vdtheta(exogct)) then
								dfpolind(exogct,endogct)=0
							else
								dfpolind(exogct,endogct)=1
							end if
							if (bct==1) dfpolind(exogct,endogct)=0
							if (dfpolind(exogct,endogct)==0) then
								Vtheta(exogct,endogct)=bestvalue
							else
								Vtheta(exogct,endogct)=Vdtheta(exogct)
							end if
							lazy_gap(exogct,endogct) = &
								lazy_gap(exogct,endogct)-lazy_range(exogct)
						end do
						envelope_candidates = envelope_candidates+bnum
						cycle
					end if
					posbinval(:) = -huge(1.0d0)
					negbinval(:) = -huge(1.0d0)
					posbincount(:) = 0
					negbincount(:) = 0
					poshead(:) = 0
					postail(:) = 0
					neghead(:) = 0
					negtail(:) = 0
					posnext(:) = 0
					negnext(:) = 0

					do kbpolct=bnum+1,endognum
						kprimect = endogind(kbpolct,1)
						bprimeval = endogval(kbpolct,2)
						qval = qmatT(kbpolct,exogct)
						actionintercept = (1.0d0-tau)*&
							(ynetcache(Act,zct,kct)-phi) &
							+ invprice*(1.0d0-delta)*kval + tau*delta*kval &
							+ actionkbase(kct,kprimect) + qval*bprimeval

						! Raw dividends fall monotonically over the ordered inherited-debt
						! grid. Locate the last nonnegative node in O(log bnum), using
						! the same floating-point expression as the former linear scan.
						splitct = 0
						splitlo = 1
						splithi = bnum
						do while (splitlo<=splithi)
							splitmid = (splitlo+splithi)/2
							rawdividend = actionintercept + dividendslope*b0(splitmid)
							if (rawdividend>=0.0d0) then
								splitct = splitmid
								splitlo = splitmid+1
							else
								splithi = splitmid-1
							end if
						end do

						posconstant = actionintercept + EVthetaT(kbpolct,exogct)
						negconstant = (1.0d0+etad)*actionintercept-etaf &
							+ EVthetaT(kbpolct,exogct)
						actionsplit(kbpolct) = splitct
						posactionval(kbpolct) = posconstant
						negactionval(kbpolct) = negconstant
						if (posconstant>posbinval(splitct)) then
							posbinval(splitct) = posconstant
						end if
						if (negconstant>negbinval(splitct)) then
							negbinval(splitct) = negconstant
						end if
					end do

					! Retain the same 1e-12 numerical tie guard used by the current
					! interval-guard authority, now within each exact-arithmetic bin.
					do kbpolct=bnum+1,endognum
						splitct = actionsplit(kbpolct)
						if (posactionval(kbpolct)>=posbinval(splitct)-tie_guard) then
							if (posbincount(splitct)==0) then
								poshead(splitct) = kbpolct
							else
								posnext(postail(splitct)) = kbpolct
							end if
							postail(splitct) = kbpolct
							posbincount(splitct) = posbincount(splitct)+1
						end if
						if (negactionval(kbpolct)>=negbinval(splitct)-tie_guard) then
							if (negbincount(splitct)==0) then
								neghead(splitct) = kbpolct
							else
								negnext(negtail(splitct)) = kbpolct
							end if
							negtail(splitct) = kbpolct
							negbincount(splitct) = negbincount(splitct)+1
						end if
					end do

					do bct=1,bnum
						endogct = bnum*(kct-1)+bct
						bval = b0(bct)
						currentdebtshift = dividendslope*bval
						bestapprox = -huge(1.0d0)
						do splitct=bct,bnum
							if (posbincount(splitct)>0) then
								regimeapprox = posbinval(splitct)+currentdebtshift
								bestapprox = max(bestapprox,regimeapprox)
							end if
						end do
						do splitct=0,bct-1
							if (negbincount(splitct)>0) then
								regimeapprox = negbinval(splitct) &
									+ (1.0d0+etad)*currentdebtshift
								bestapprox = max(bestapprox,regimeapprox)
							end if
						end do

						ntop = 0
						! Eligible split bins are disjoint between the positive and
						! issuance regimes. Append only global near-winner bin tie sets;
						! never rescan all 870 actions for each of the 30 debt states.
						do splitct=bct,bnum
							regimeapprox = posbinval(splitct)+currentdebtshift
							if (posbincount(splitct)>0 .and. &
								regimeapprox>=bestapprox-tie_guard) then
								candidatect = poshead(splitct)
								do while (candidatect/=0)
									ntop = ntop+1
									topind(ntop) = candidatect
									candidatect = posnext(candidatect)
								end do
							end if
						end do
						do splitct=0,bct-1
							regimeapprox = negbinval(splitct) &
								+ (1.0d0+etad)*currentdebtshift
							if (negbincount(splitct)>0 .and. &
								regimeapprox>=bestapprox-tie_guard) then
								candidatect = neghead(splitct)
								do while (candidatect/=0)
									ntop = ntop+1
									topind(ntop) = candidatect
									candidatect = negnext(candidatect)
								end do
							end if
						end do

						do topct=1,ntop
							kbpolct = topind(topct)
							kprimect = endogind(kbpolct,1)
							bprimeval = endogval(kbpolct,2)
							qval = qmat(exogct,kbpolct)
							dval = dividend(A0(Act),z0(zct),kval,k0(kprimect),&
								qval,bprimeval,bval)
							sparsevalue(topct) = dval
							sparsedividend(topct) = dval
						end do
						envelope_candidates = envelope_candidates+ntop

						do Aprimect=1,Anum
							do zprimect=1,znum
								exogprimect = znum*(Aprimect-1)+zprimect
								do topct=1,ntop
									kbpolct = topind(topct)
									if (doRATfirms==0) then
										sparsevalue(topct) = sparsevalue(topct) + beta * &
											pr_mat_exog_theta(exogct,exogprimect) * &
											Vthetaold(exogprimect,kbpolct)
									else if (doRATfirms==1) then
										sparsevalue(topct) = sparsevalue(topct) + beta * &
											pr_mat_exog(exogct,exogprimect) * &
											Vthetaold(exogprimect,kbpolct)
									end if
								end do
							end do
						end do

						sparse_nonfinite = .false.
						do topct=1,ntop
							if (.not.ieee_is_finite(sparsevalue(topct))) &
								sparse_nonfinite = .true.
						end do
						if (sparse_nonfinite) then
							RHSndvec(:) = -1000000.0d0
							RHSdivvec(:) = -1000000.0d0
							do topct=1,ntop
								kbpolct = topind(topct)
								RHSndvec(kbpolct) = sparsevalue(topct)
								RHSdivvec(kbpolct) = sparsedividend(topct)
							end do
							bestaction = maxloc(RHSndvec,1)
							bestvalue = RHSndvec(bestaction)
							bestdividend = RHSdivvec(bestaction)
						else
							bestaction = topind(1)
							bestvalue = sparsevalue(1)
							bestdividend = sparsedividend(1)
							do topct=2,ntop
								kbpolct = topind(topct)
								if (sparsevalue(topct)>bestvalue .or. &
									(sparsevalue(topct)==bestvalue .and. kbpolct<bestaction)) then
									bestaction = kbpolct
									bestvalue = sparsevalue(topct)
									bestdividend = sparsedividend(topct)
								end if
							end do
						end if
						lazy_secondbest = -huge(1.0d0)
						do topct=1,ntop
							if (topind(topct)/=bestaction) &
								lazy_secondbest = max(lazy_secondbest,sparsevalue(topct))
						end do
						if (ntop==1) then
							lazy_margin = 0.5d0*tie_guard
						else
							lazy_margin = min(0.5d0*tie_guard,&
								max(0.0d0,bestvalue-lazy_secondbest))
						end if
						lazy_gap(exogct,endogct) = lazy_margin
						kbpolct = bestaction
						kbpolind(exogct,endogct) = kbpolct
						Vndtheta(exogct,endogct) = bestvalue
						if (bestdividend>=0.0d0) then
							divposmat(exogct,endogct) = 1.0d0
						else
							divposmat(exogct,endogct) = 0.0d0
						end if

						if (Vndtheta(exogct,endogct)>=Vdtheta(exogct)) then
							dfpolind(exogct,endogct)=0
						else
							dfpolind(exogct,endogct)=1
						end if
						if (bct==1) dfpolind(exogct,endogct)=0
						if (dfpolind(exogct,endogct)==0) then
							Vtheta(exogct,endogct)=Vndtheta(exogct,endogct)
						else
							Vtheta(exogct,endogct)=Vdtheta(exogct)
						end if
					end do
				end do
			end do
			!$omp end parallel do
			lazy_prev_ev(:,:) = EVthetaT(:,:)
			lazy_valid = .true.
			write(*,'(A,I0,A,I0,A,F10.6)') 'RCC_LAZY_BLOCKS=',lazy_blocks,&
				'TOTAL=',lazy_total_blocks,' SHARE=',&
				dble(lazy_blocks)/dble(max(1_8,lazy_total_blocks))
			write(*,'(A,I0,A,I0,A,F10.6)') 'RCC_LAZY_CERTIFIED_STATES=',&
				lazy_certified_states,' TOTAL=',lazy_total_states,' SHARE=',&
				dble(lazy_certified_states)/dble(max(1_8,lazy_total_states))
			write(*,'(A,I0,A,F10.4)') 'RCC_ENVELOPE_CANDIDATES=',envelope_candidates,&
				' MEAN=',dble(envelope_candidates)/dble(exognum*endognum)

			!now, compute & report the implied error in policies & value functions
			polerr = maxval(abs(kbpolind-kbpolindold))
			vferr = maxval(abs(Vtheta-Vthetaold))
			write(*,*) "Pol Ct = ",polct," Pol Err = ",polerr," V Err = ",vferr

			!if policies converged, then exit
			if (polerr<1.0) then
      !if (vferr<vftol) then
				solveFLAG = 0
				exit
			end if

			!if not, then update and continue
			kbpolindold = kbpolind

      minshift = (beta/(1.0-beta))*minval(Vtheta-Vthetaold);
      maxshift = (beta/(1.0-beta))*maxval(Vtheta-Vthetaold);
      !write(*,*) "min = ",minshift," max = ",maxshift

      Vthetaold = Vtheta

		end do !polct

    call distribution()

		!check for convergence of policies with outer loop policies
		dferr = maxval(abs(dfpolind - dfpolindold))
    !dfmeanerr = sum(abs(dfpolind - dfpolindold))/dble(exognum*endognum)
    dfmeanerr = sum(dist*abs(dfpolind-dfpolindold))
		write(*,*) "Implied max dferr = ",dferr," mean dferr = ",dfmeanerr

		!if (dferr<1.0) then
    !if ((dfmeanerr<poltol).and.(dfct>3)) then
    if ((dfmeanerr<poltol).and.(dfct>1)) then
			solveFLAG = 0
			exit
		end if

		dfpolindold = dfpolind


	end do !dfct

!or do yo just read in the policies and values?
else if (doVFI==0) then

	write(*,*) "##### Not doing policy iteration "
	write(*,*) "Reading in policies + values instead of computing them."

	!value function
	open(13,file="Vtheta.txt")
	do exogct=1,exognum
	do endogct=1,endognum
		read(13,*) Vtheta(exogct,endogct)
	end do !endogct
	end do !exogct
	close(13)

	!policy function
	open(13,file="kbpolind.txt")
	do exogct=1,exognum
	do endogct=1,endognum
		read(13,*) kbpolind(exogct,endogct)
	end do !endogct
	end do !exogct
	close(13)

	!default function
	open(13,file="dfpolind.txt")
	do exogct=1,exognum
	do endogct=1,endognum
		read(13,*) dfpolind(exogct,endogct)
	end do !endogct
	end do !exogct
	close(13)

	!bond price function
	open(13,file="qmat.txt")
	do exogct=1,exognum
	do endogct=1,endognum
		read(13,*) qmat(exogct,endogct)
	end do !endogct
	end do !exogct
	close(13)

	write(*,*) "Done with read-in process."
	write(*,*) " "



end if

if (doVFI==1 .and. lazy_total_blocks>0) then
	write(*,'(A,I0,A,I0,A,F12.8)') 'RCC_LAZY_SUMMARY_SKIPPED=',lazy_blocks,&
		' TOTAL=',lazy_total_blocks,' SHARE=',&
		dble(lazy_blocks)/dble(max(1_8,lazy_total_blocks))
	write(*,'(A,I0,A,I0,A,F12.8)') 'RCC_LAZY_SUMMARY_CERTIFIED_STATES=',&
		lazy_certified_states,' TOTAL=',lazy_total_states,' SHARE=',&
		dble(lazy_certified_states)/dble(max(1_8,lazy_total_states))
end if

if (doVFI==1 .and. active_updates>0) then
	write(*,'(A,I0,A,F12.6,A,I0,A,I0)') 'RCC_ACTIVE_PUSH_UPDATES=',active_updates,&
		'MEAN_ACTIONS=',dble(active_total)/dble(active_updates),&
		' MIN_ACTIONS=',active_min,' MAX_ACTIONS=',active_max
end if

deallocate(lazy_prev_ev,lazy_gap,lazy_range)
deallocate(EVtheta,EVthetaT,qmatT,currentdiv,ynetcache,accache,currentbase,&
	actionkbase,debtcost,actioncost,actionb,actioninvalid)

!GMM diagnostics
solveFLAG=0
if (dfct==(dfmaxit-1)) then
  solveFLAG=1
end if

!determine the profit matrix
!NOTE: tried parallelization here, and no gain due to overhead
do exogct=1,exognum
do endogct=1,endognum

  !extract states
  Amin1val = exogval(exogct,1)
  Aval = exogval(exogct,2)
  zmin1val = exogval(exogct,3)
  zval = exogval(exogct,4)
  Act = exogind(exogct,2)
  zct = exogind(exogct,4)
  kval = endogval(endogct,1)
  bval = endogval(endogct,2)

  dfct = dfpolind(exogct,endogct)

  if (dfct==0) then

    kbpolct = kbpolind(exogct,endogct)

  else if (dfct==1) then

    !dumkvec = (k0-(1.0-delta)*endogval(endogct,1)*gamma)**2.0
    kprimect = 1
    dumbvec = b0**2.0
    bprimect = minloc(dumbvec,1)
    kbpolct = bnum*(kprimect-1) + bprimect

  end if

  kprimeval = endogval(kbpolct,1)
  bprimeval = endogval(kbpolct,2)
  yval = output(Aval,zval,kval) * (1.0-dble(dfct)*(1.0-gamma))
  ival = kprimeval - (1.0-delta)*kval
  !qval = qmat(exogct,endogct)
  qval = qmat(exogct,kbpolct)
  ACval = AC(kval,kprimeval)

  pimat(exogct,endogct) = (1.0-tau)*(outputnetlabor(Aval,zval,kval)  &
		- ACval - phi) - delta*kval + tau*(r*bval + delta*kval)

end do !endogct
end do !exogct

end subroutine politeration


subroutine discretize()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,zct,zprimect,zmin1ct,Amin1ct,Act,Aprimect,exogct,exogprimect,&
	bct,kct
double precision :: zmin1val,zval,zprimeval,newsval,gridinc,nextmean,&
	Amin1val,Aval,Aprimeval,bval,kval,ergdistAnew(Anum),ergdistAold(Anum),k0intensive(knum-1),&
  b0intensive(bnum-1)

write(*,*) "##### Discretizing the model "

!!!!!!!discretize the capital grid with log-linear spacing
call linspace(k0intensive,log(kmin),log(kmax),knum-1)
k0intensive = exp(k0intensive)
!call linspace(k0intensive,kmin,kmax,knum-1)
k0(2:knum) = k0intensive
k0(1) = 0.0

!!!!!!!discretize the debt grid with linear spacing
call linspace(b0intensive,log(bmin),log(bmax),bnum-1)
b0intensive = exp(b0intensive)
!call linspace(b0intensive,bmin,bmax,bnum-1)
b0(2:bnum) = b0intensive
b0(1) = 0.0


!!!!!!!now, discretize the true process for z
call idio_tauchen(znum,rhoz,sigmaz,pr_mat_z,z0,nsigmaz)

!!!!!!!now, discretize the true process for A
call idio_tauchen(Anum,rhoA,sigmaA,pr_mat_A,A0,nsigmaA)

!compute the implied ergodic distribution of A
ergdistAold(:) = 1.0/dble(Anum)
do ct=1,distmaxit

  ergdistAnew(:) = 0.0

  do Act=1,Anum

    do Aprimect=1,Anum
      ergdistAnew(Aprimect) = ergdistAnew(Aprimect) + pr_mat_A(Act,Aprimect)*&
        ergdistAold(Act)
    end do !Aprimect

  end do !Act

  ergdistAnew = ergdistAnew/sum(ergdistAnew)
  if (maxval(abs(ergdistAnew - ergdistAold))<disttol) then
    exit
  end if

  ergdistAold = ergdistAnew

end do !ct

ergdistA = ergdistAold

!!!!!!!discretize the spread grid with linear spacing
call linspace(spread0,spreadmin,spreadmax,spreadnum)

!!!!!!!discretize the diagnostic expectations process for idio prod

!first, convert the grid back to logs from levels
z0 = log(z0)
gridinc = z0(2)-z0(1)

!then, do the construction of the prob matrix
pr_mat_z_theta(:,:,:) = 0.0 !(z_{-1},z,z')
do zmin1ct=1,znum

	do zct=1,znum

		!determine the mean of next period's diagnostic dist
		zmin1val = z0(zmin1ct)
		zval = z0(zct)
		newsval = rhoz * (zval - rhoz * zmin1val)
		nextmean = rhoz * zval + theta * newsval

		!insert transition matrix middle rows
		do zprimect=2,(znum-1)

			pr_mat_z_theta(zmin1ct,zct,zprimect) = &
				normcdf(z0(zprimect)+gridinc/2.0,nextmean,sigmaz) - &
				normcdf(z0(zprimect)-gridinc/2.0,nextmean,sigmaz)

		end do !zprimect

	    !first interval and last interval take the rest of the weight
	    pr_mat_z_theta(zmin1ct,zct,1) = normcdf(z0(1)+gridinc/2.0,nextmean,sigmaz)
	    pr_mat_z_theta(zmin1ct,zct,znum) = 1.0 - normcdf(z0(znum)-gridinc/2.0,nextmean,sigmaz)


end do !zct
end do !zmin1ct

!convert the grid back to levels rather than logs
z0 = exp(z0)

!!!!!!!discretize the diagnostic expectations process for agg prod

!first, convert the grid back to logs from levels
A0 = log(A0)
gridinc = A0(2)-A0(1)

!then, do the construction of the prob matrix
pr_mat_A_theta(:,:,:) = 0.0 !(A_{-1},A,A')
do Amin1ct=1,Anum


	do Act=1,Anum

		!determine the mean of next period's diagnostic dist
		Amin1val = A0(Amin1ct)
		Aval = A0(Act)
		newsval = rhoA * (Aval - rhoA * Amin1val)
		nextmean = rhoA * Aval + theta * newsval

		!insert transition matrix middle rows
		do Aprimect=2,(Anum-1)

			pr_mat_A_theta(Amin1ct,Act,Aprimect) = &
				normcdf(A0(Aprimect)+gridinc/2.0,nextmean,sigmaA) - &
				normcdf(A0(Aprimect)-gridinc/2.0,nextmean,sigmaA)

		end do !Aprimect

	    !first interval and last interval take the rest of the weight
	    pr_mat_A_theta(Amin1ct,Act,1) = normcdf(A0(1)+gridinc/2.0,nextmean,sigmaA)
	    pr_mat_A_theta(Amin1ct,Act,Anum) = 1.0 - normcdf(A0(Anum)-gridinc/2.0,nextmean,sigmaA)


end do !zct
end do !zmin1ct

!convert the grid back to levels rather than logs
A0 = exp(A0)

write(*,*) " "


!!!!!!!now, do the joint discretization into a single exogenous variable

!Under theta=0, lagged A and z labels are redundant.  Retain the
!four-column interface used by the released routines, but duplicate current
!indices into the lag columns and allocate only the 5 x 5 sufficient states.
ct = 0
do Act=1,Anum
  do zct=1,znum
    ct = ct + 1
    exogind(ct,1) = Act
    exogind(ct,2) = Act
    exogind(ct,3) = zct
    exogind(ct,4) = zct
    exogval(ct,1) = A0(Act)
    exogval(ct,2) = A0(Act)
    exogval(ct,3) = z0(zct)
    exogval(ct,4) = z0(zct)
  end do
end do

!now, create the joint transition matrix
pr_mat_exog(:,:) = 0.0
pr_mat_exog_theta(:,:) = 0.0

do exogct=1,exognum

	!extract current states
	Amin1ct = exogind(exogct,1)
	Act = exogind(exogct,2)
	zmin1ct = exogind(exogct,3)
	zct = exogind(exogct,4)

	!now, loop over next period realizations Aprimect & zprimect
	do Aprimect=1,Anum
		do zprimect=1,znum
			exogprimect = znum*(Aprimect-1) + zprimect

			pr_mat_exog(exogct,exogprimect) = pr_mat_A(Act,Aprimect)*pr_mat_z(zct,zprimect)

			pr_mat_exog_theta(exogct,exogprimect) = pr_mat_A_theta(Amin1ct,Act,Aprimect)&
				*pr_mat_z_theta(zmin1ct,zct,zprimect)

		end do !zprimect
	end do !Aprimect

	!correct for rounding error
	pr_mat_exog(exogct,:) = 	pr_mat_exog(exogct,:)  / sum(	pr_mat_exog(exogct,:)  )
	pr_mat_exog_theta(exogct,:) = 	pr_mat_exog_theta(exogct,:)  / sum(	pr_mat_exog_theta(exogct,:)  )

end do !exogct

!now, do the indexing for the endogenous states
endogval(:,:) = 0.0
endogind(:,:) = 0

ct = 0
do kct=1,knum
	do bct=1,bnum
		ct = ct + 1
		kval = k0(kct)
		bval = b0(bct)

		endogval(ct,1) = kval
		endogval(ct,2) = bval

		endogind(ct,1) = kct
		endogind(ct,2) = bct

	end do !bct
end do !kct


end subroutine discretize

subroutine writefiles()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,zmin1ct,zct,kct,zprimect,firmct,t,Amin1ct,Act,Aprimect,ct2,exogct,Kbarct,bct,&
	endogct,simct,shockct

write(*,*) "##### Writing output files "

!constants
open(13,file="constants.txt")
write(13,*) (/dble(znum),dble(Anum),dble(knum),dble(bnum),dble(exognum),dble(endognum),&
	dble(Tsim),dble(Terg),dble(Tirf),dble(Nirf),dble(TstartIRF),&
  dble(Tnumirf),dble(doIRF),dble(Nfirm),dble(spreadnum)/)
close(13)

!z0 grid
open(13,file="z0.txt")
do ct=1,znum
	write(13,*) z0(ct)
end do !ct
close(13)

!z0 grid
open(13,file="spread0.txt")
do ct=1,spreadnum
	write(13,*) spread0(ct)
end do !ct
close(13)

!true pr_mat_z transitions
open(13,file="pr_mat_z.txt")
do ct=1,znum
	write(13,*) pr_mat_z(ct,:)
end do !ct
close(13)

!diagnostic expectations transitions
open(13,file="pr_mat_z_theta.txt")
do zmin1ct=1,znum
do zct=1,znum
do zprimect=1,znum
	write(13,*) pr_mat_z_theta(zmin1ct,zct,zprimect)
end do !ct3
end do !ct2
end do !ct1
close(13)

!A0 grid
open(13,file="A0.txt")
do ct=1,Anum
	write(13,*) A0(ct)
end do !ct
close(13)

!true pr_mat_A transitions
open(13,file="pr_mat_A.txt")
do ct=1,Anum
	write(13,*) pr_mat_z(ct,:)
end do !ct
close(13)

!diagnostic expectations transitions
open(13,file="pr_mat_A_theta.txt")
do Amin1ct=1,Anum
do Act=1,Anum
do Aprimect=1,Anum
	write(13,*) pr_mat_A_theta(Amin1ct,Act,Aprimect)
end do !ct3
end do !ct2
end do !ct1
close(13)

!exogenous indexes
open(13,file="exogind.txt")
do ct=1,exognum
	write(13,*) exogind(ct,:)
end do !ct
close(13)

!exogenous values
open(13,file="exogval.txt")
do ct=1,exognum
	write(13,*) exogval(ct,:)
end do !ct
close(13)

!k0 grid
open(13,file="k0.txt")
do ct=1,knum
	write(13,*) k0(ct)
end do !ct
close(13)

!kdense0 grid
open(13,file="b0.txt")
do ct=1,bnum
	write(13,*) b0(ct)
end do !ct
close(13)

!endogenous values
open(13,file="endogval.txt")
do ct=1,endognum
	write(13,*) endogval(ct,:)
end do !ct
close(13)

!endogenous indexes
open(13,file="endogind.txt")
do ct=1,endognum
	write(13,*) endogind(ct,:)
end do !ct
close(13)


if (doVFI==1) then
!value function
open(13,file="Vtheta.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) Vtheta(exogct,endogct)
end do !endogct
end do !exogct
close(13)

!value no default
open(13,file="Vndtheta.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) Vndtheta(exogct,endogct)
end do !endogct
end do !exogct
close(13)

!value default
open(13,file="Vdtheta.txt")
do exogct=1,exognum
	write(13,*) Vdtheta(exogct)
end do !exogct
close(13)

!policy function
open(13,file="kbpolind.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) kbpolind(exogct,endogct)
end do !endogct
end do !exogct
close(13)

!default function
open(13,file="dfpolind.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) dfpolind(exogct,endogct)
end do !endogct
end do !exogct
close(13)

!bond price function
open(13,file="qmat.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) qmat(exogct,endogct)
end do !endogct
end do !exogct
close(13)

end if

!distribution
open(13,file="dist.txt")
do exogct=1,exognum
do endogct=1,endognum
	write(13,*) dist(exogct,endogct)
end do !endogct
end do !exogct
close(13)


!distribution k
open(13,file="distk.txt")
do ct=1,knum
	write(13,*) distk(ct)
end do !ct
close(13)

!distribution b
open(13,file="distb.txt")
do ct=1,bnum
	write(13,*) distb(ct)
end do !ct
close(13)

!Asim
open(13,file="Asim.txt")
do t=1,Tsim
write(13,*) Asim(t)
end do !t
close(13)

!Asimpos
open(13,file="Asimpos.txt")
do t=1,Tsim
write(13,*) Asimpos(t)
end do !t
close(13)

!Ysim
open(13,file="Ysim.txt")
do t=1,Tsim
write(13,*) Ysim(t)
end do !t
close(13)

!Ksim
open(13,file="Ksim.txt")
do t=1,Tsim
write(13,*) Ksim(t)
end do !t
close(13)

!Isim
open(13,file="Isim.txt")
do t=1,Tsim
write(13,*) Isim(t)
end do !t
close(13)

!qsim
open(13,file="qsim.txt")
do t=1,Tsim
write(13,*) qsim(t)
end do !t
close(13)

!Bsim
open(13,file="Bsim.txt")
do t=1,Tsim
write(13,*) Bsim(t)
end do !t
close(13)

!dfsim
open(13,file="dfsim.txt")
do t=1,Tsim
write(13,*) dfsim(t)
end do !t
close(13)

!BtoKsim
open(13,file="BtoKsim.txt")
do t=1,Tsim
write(13,*) BtoKsim(t)
end do !t
close(13)

!spreadsim
open(13,file="spreadsim.txt")
do t=1,Tsim
write(13,*) spreadsim(t)
end do !t
close(13)

!spreadBAAsim
open(13,file="spreadBAAsim.txt")
do t=1,Tsim
write(13,*) spreadBAAsim(t)
end do !t
close(13)

!spreadsim
open(13,file="retbondsim.txt")
do t=1,Tsim
write(13,*) retbondsim(t)
end do !t
close(13)

!spreadtracksim
open(13,file="spreadtracksim.txt")
do ct=1,spreadnum
do t=1,Tsim
write(13,*) spreadtracksim(ct,t)
end do !t
end do !ct
close(13)

!Nsim
open(13,file="Nsim.txt")
do t=1,Tsim
write(13,*) Nsim(t)
end do !t
close(13)

!sdMPKsim
open(13,file="sdMPKsim.txt")
do t=1,Tsim
write(13,*) sdMPKsim(t)
end do !t
close(13)

!issjunksim
open(13,file="issjunksim.txt")
do t=1,Tsim
write(13,*) issjunksim(t)
end do !t
close(13)


!dfYsim
open(13,file="dfYsim.txt")
do t=1,Tsim
write(13,*) dfYsim(t)
end do !t
close(13)

!dfIsim
open(13,file="dfIsim.txt")
do t=1,Tsim
write(13,*) dfIsim(t)
end do !t
close(13)

!dfBsim
open(13,file="dfBsim.txt")
do t=1,Tsim
write(13,*) dfBsim(t)
end do !t
close(13)

!FEsim
open(13,file="FEsim.txt")
do t=1,Tsim
write(13,*) FEsim(t)
end do !t
close(13)

!Pisim
open(13,file="Pisim.txt")
do t=1,Tsim
write(13,*) Pisim(t)
end do !t
close(13)

!Pisim
open(13,file="Pifsim.txt")
do t=1,Tsim
write(13,*) Pifsim(t)
end do !t
close(13)

!irateheterosim
open(13,file="irateheterosim.txt")
do ct=1,6
do t=1,Tsim
write(13,*) irateheterosim(t,ct)
end do !t
end do !ct
close(13)

!Dbheterosim
open(13,file="Dbheterosim.txt")
do ct=1,6
do t=1,Tsim
write(13,*) Dbheterosim(t,ct)
end do !t
end do !ct
close(13)

!zfirmpos
open(13,file="zfirmpos.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) zfirmpos(t,firmct)
end do !t
end do !firmct
close(13)

!endogfirmpos
open(13,file="endogfirmpos.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) endogfirmpos(t,firmct)
end do !t
end do !firmct
close(13)

!dfirmpos
open(13,file="dfirmpos.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) dfirmpos(t,firmct)
end do !t
end do !firmct
close(13)

!kfirmpos
open(13,file="kfirmpos.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) kfirmpos(t,firmct)
end do !t
end do !firmct
close(13)

!bfirmpos
open(13,file="bfirmpos.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) bfirmpos(t,firmct)
end do !t
end do !firmct
close(13)

!pifirm
open(13,file="pifirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) pifirm(t,firmct)
end do !t
end do !firmct
close(13)

!pifcstfirm
open(13,file="pifcstfirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) pifcstfirm(t,firmct)
end do !t
end do !firmct
close(13)

!ifirm
open(13,file="ifirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) ifirm(t,firmct)
end do !t
end do !firmct
close(13)

!kfirm
open(13,file="kfirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) kfirm(t,firmct)
end do !t
end do !firmct
close(13)

!bfirm
open(13,file="bfirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) bfirm(t,firmct)
end do !t
end do !firmct
close(13)

!fefirm
open(13,file="fefirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) fefirm(t,firmct)
end do !t
end do !firmct
close(13)

!qfirm
open(13,file="qfirm.txt")
do firmct=1,Nfirm
do t=1,Tsim
write(13,*) qfirm(t,firmct)
end do !t
end do !firmct
close(13)

!do the IRF simulation output
if (doIRF==1) then

!AsimposIRF
open(13,file="AsimposIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) AsimposIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!AsimIRF
open(13,file="AsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) AsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!YsimIRF
open(13,file="YsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) YsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!KsimIRF
open(13,file="KsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) KsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!IsimIRF
open(13,file="IsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) IsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!qsimIRF
open(13,file="qsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) qsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!BsimIRF
open(13,file="BsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) BsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!dfsimIRF
open(13,file="dfsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) dfsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!BtoKsimIRF
open(13,file="BtoKsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) BtoKsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!spreadsimIRF
open(13,file="spreadsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) spreadsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!FEsimIRF
open(13,file="FEsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) FEsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!PisimIRF
open(13,file="PisimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) PisimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

!PifsimIRF
open(13,file="PifsimIRF.txt")
do shockct=1,2
do simct=1,Nirf
do t=1,Tirf
write(13,*) PifsimIRF(t,simct,shockct)
end do !t
end do !simct
end do !shockct
close(13)

end if !doIRF

if (doMATCH==1) then

open(13,file="Amatch.txt")
do t=1,Tmatch
write(13,*) Amatch(t)
end do !t
close(13)

open(13,file="Amatchpos.txt")
do t=1,Tmatch
write(13,*) Amatchpos(t)
end do !t
close(13)

open(13,file="pAmatch.txt")
do t=1,Tmatch
write(13,*) pAmatch(t)
end do !t
close(13)

open(13,file="Ymatch.txt")
do t=1,Tmatch
write(13,*) Ymatch(t)
end do !t
close(13)

open(13,file="Imatch.txt")
do t=1,Tmatch
write(13,*) Imatch(t)
end do !t
close(13)

open(13,file="dImatch.txt")
do t=1,Tmatch
write(13,*) dImatch(t)
end do !t
close(13)

open(13,file="qmatch.txt")
do t=1,Tmatch
write(13,*) qmatch(t)
end do !t
close(13)

open(13,file="dfmatch.txt")
do t=1,Tmatch
write(13,*) dfmatch(t)
end do !t
close(13)

open(13,file="Bmatch.txt")
do t=1,Tmatch
write(13,*) Bmatch(t)
end do !t
close(13)

open(13,file="BtoKmatch.txt")
do t=1,Tmatch
write(13,*) BtoKmatch(t)
end do !t
close(13)

open(13,file="Kmatch.txt")
do t=1,Tmatch
write(13,*) Kmatch(t)
end do !t
close(13)

open(13,file="Nmatch.txt")
do t=1,Tmatch
write(13,*) Nmatch(t)
end do !t
close(13)

open(13,file="sdMPKmatch.txt")
do t=1,Tmatch
write(13,*) sdMPKmatch(t)
end do !t
close(13)

open(13,file="spreadmatch.txt")
do t=1,Tmatch
write(13,*) spreadmatch(t)
end do !t
close(13)

open(13,file="spreadBAAmatch.txt")
do t=1,Tmatch
write(13,*) spreadBAAmatch(t)
end do !t
close(13)

open(13,file="dspreadBAAmatch.txt")
do t=1,Tmatch
write(13,*) dspreadBAAmatch(t)
end do !t
close(13)

open(13,file="dspreadmatch.txt")
do t=1,Tmatch
write(13,*) dspreadmatch(t)
end do !t
close(13)

end if !doMATCH


write(*,*) " "



end subroutine writefiles

double precision function AC(kval,kprimeval)
use params
use base_lib
use omp_lib
implicit none

double precision :: kval,kprimeval,ival

ival = kprimeval - (1.0-delta)*kval

if (kval > 0.0) then
AC = (etak/2.0)
AC = AC * (ival/kval)**2.0
AC = AC * kval
else
  AC = 0.0
end if

end function AC


double precision function labor(Aval,zval,kval)
use params
use base_lib
use omp_lib
implicit none

double precision :: Aval,zval,kval

labor = ((nu/w)**(1.0/(1.0-nu)))* &
	((Aval*zval)**(1.0/(1.0-nu)))* (kval**(alpha/(1.0-nu)))

end function labor

double precision function output(Aval,zval,kval)
use params
use base_lib
use omp_lib
implicit none

double precision :: Aval,zval,kval

output = ((nu/w)**(nu/(1.0-nu)))* ((zval*Aval)**(1.0/(1.0-nu)))* &
	(kval**(alpha/(1.0-nu)))

end function output

double precision function outputnetlabor(Aval,zval,kval)
use params
use base_lib
use omp_lib
implicit none

double precision :: Aval,zval,kval

outputnetlabor = (1.0-nu)*((nu/w)**(nu/(1.0-nu)))* &
	((Aval*zval)**(1.0/(1.0-nu)))* (kval**(alpha/(1.0-nu)))

end function outputnetlabor

double precision function dividend(Aval,zval,kval,kprimeval,&
	qval,bprimeval,bval)

double precision :: Aval,zval,kval,kprimeval,ival,ACval,qval,&
	bprimeval,bval

ival = kprimeval - (1.0-delta)*kval
ACval = AC(kval,kprimeval)

	dividend = (1.0-tau)*(outputnetlabor(Aval,zval,kval)  &
		- ACval - phi) - invprice*ival + tau*(r*bval + delta*kval) + &
		qval * bprimeval - bval

  if (bprimeval>bval) then
    dividend = dividend - (etab/2.0)*(bprimeval-bval)**dble(2.0)
  end if

  if (dividend<0.0) then
    dividend = dividend - etaf - etad*abs(dividend)
  end if

end function dividend

double precision function recovery(kval,bval,Aval,zval)

double precision :: kval,bval,Aval,zval,ynetwval

ynetwval = outputnetlabor(Aval,zval,kval)

if (bval<=0.0) then

  recovery = 1.0

else
   recovery =  (1-tau)*gamma*((1.0-delta)*kval)/bval

end if

if (recovery>1.0) then

	recovery=1.0

end if

if (recovery<=0.0) then

	recovery=0.0

end if


end function recovery

end module solve
