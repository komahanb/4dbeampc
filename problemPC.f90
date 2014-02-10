program problemPCBeam
  use dimopt
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt

  implicit none
  !
  !     include the Ipopt return codes
  !
  include 'IpReturnCodes.inc'
  include 'mpif.h'
  !
  !     Size of the problem (number of variables and equality constraints)
  !
  integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
  parameter  (N = 4, M = 3, NELE_JAC = 12, NELE_HESS = 10)
  parameter  (IDX_STY = 1 )
  !
  !     Space for multipliers and constraints
  !
  double precision LAM(M)
  double precision G(M)
  !
  !     Vector of variables
  !
  double precision X(N)
  !
  !     Vector of lower and upper bounds
  !
  double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
  double precision G_L(M), G_U(M)
  !
  !     Private data for evaluation routines
  !     This could be used to pass double precision and integer arrays untouched
  !     to the evaluation subroutines EVAL_*
  !
  double precision DAT(2000)
  integer IDAT(2000)
  !
  !     Place for storing the Ipopt Problem Handle
  !
  integer*8 IPROBLEM
  integer*8 IPCREATE
  !
  integer IERR
  integer IPSOLVE, IPADDSTROPTION
  integer IPADDNUMOPTION, IPADDINTOPTION
  integer IPOPENOUTPUTFILE
  !
  double precision F,Fs,sigmax(N)
  integer i,kprob

  double precision  infbound
  parameter        (infbound = 1.d+20)
  !
  !     The following are the Fortran routines for computing the model
  !     functions and their derivatives - their code can be found further
  !     down in this file.
  !
  external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS, ITER_CB


  call MPI_START

  !**********************************************************************


!  pi=4.0*atan(1.0) ! constant for later use (visible globally)

  !======================================
  !(1)    Set initial point and bounds:
  !=======================================

  ! Evals counter initialization

  fcnt=0
  fgcnt=0
  fghcnt=0

  !Breadth and depth of the bar

  do i=1,N-2
     X(i)   = 300.0  
     X_L(i) = 100.0 
     X_U(i) = 600.0
  end do

  ! Bending moment (Aleatory)

  X(3)  = 40.0e6 !N.m
  X_L(3)= 40.0e6 !N.m
  X_U(3)= 40.0e6 !N.m

  ! Shear force  (Aleatory)

  X(4)  = 150.0e3 !N
  X_L(4)= 150.0e3 !N
  X_U(4)= 150.0e3 !N


  !===================================================================
  !(2)Integer Settings and store into IDAT (check for size above)
  !===================================================================

  probtype(:)=1
  kprob=0

  IDAT(1)=kprob
  IDAT(2)=0
  IDAT(3:N+2)=probtype(1:N)

  !===============================================
  !(3)     Setup std dev and store in to dat(1:N)
  !===============================================

  sigmax(1)=10.0
  sigmax(2)=10.0

  sigmax(3)=40.0d3
  sigmax(4)=1500.0

  do i=1,n
     dat(i)=sigmax(i)
  end do

  !====================
  !(4)     Constraints
  !====================

  do i=1,M
     G_L(i)=-infbound
     G_U(i)=0.d0
  end do
  
  !
  ! Equality constraint
  !
  
  G_L(M)=0.0d0
  G_U(M)=0.0d0
  
  !===========================================================
  !(5)    Other constants to be passed to the surrogate call
  !===========================================================
  
  !  pi=4.0*atan(1.0) ! constant for later use
  
  !Problem data and other constants
  
  dat(1000+1)=10.0 !Sigma_allow
  dat(1000+2)=2.0  !Tau_allow
  dat(1000+3)=1.00  !Factor of safety
  dat(1000+20)=77  !filenum for PC


  !=================================================================
  !
  !     First create a handle for the Ipopt problem (and read the options
  !     file)
  !=================================================================

  IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
  if (IPROBLEM.eq.0) then
     write(*,*) 'Error creating an Ipopt Problem handle.'
     call stop_all
  endif
  !
  !     Open an output file
  !
  if (id_proc.eq.0) open(unit=76,file='Opt.his',form='formatted',status='replace')

  if (id_proc.eq.0) open(unit=86,file='beta.his',form='formatted',status='replace')

  if (id_proc.eq.0) open(unit=37,file='dv.his',form='formatted',status='replace')


  IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
  if (IERR.ne.0 ) then
     write(*,*) 'Error opening the Ipopt output file.'
     goto 9000
  endif
  !

  !!
  !!     Set a callback function to give you control once per iteration.
  !!     You can use it if you want to generate some output, or to stop
  !!     the optimization early.
  !!
  call IPSETCALLBACK(IPROBLEM, ITER_CB)

  !
  !     Call optimization routine
  !

  if (id_proc.eq.0) then
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  else
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  end if

  IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)

  !
  !     Output:
  !
  if (id_proc.eq.0) then

     if( IERR.eq.IP_SOLVE_SUCCEEDED .or. IERR.eq.5) then
        write(*,*)
        write(*,*) 'The solution was found.'
        write(*,*)
     else
        write(*,*)
        write(*,*) 'An error occoured.'
        write(*,*) 'The error code is ',IERR
        write(*,*)
     endif

     write(*,*) 'The final value of the objective function is ',F,' m^2'
     write(*,*)
     write(*,*) 'The optimal values of X are:'
     write(*,*)
     do i = 1, N
        write(*,*) 'X  (',i,') = ',X(i)
     enddo
     write(*,*)
     write(*,*) 'The multipliers for the equality constraints are:'
     write(*,*)
     do i = 1, M
        write(*,*) 'LAM(',i,') = ',LAM(i)
     enddo
     write(*,*)
     write(*,'(a,4F13.4)') 'Weight, variance, SD, CV:',DAT(N+1),DAT(N+2),sqrt(DAT(N+2)),sqrt(DAT(N+2))/DAT(N+1)
  end if
  !
9000 continue
  !
  !     Clean up
  !
  call IPFREE(IPROBLEM)
  if (id_proc.eq.0) close(76)
  if (id_proc.eq.0) close(86)
  if (id_proc.eq.0) close(37)

  call stop_all
  !
9990 continue
  write(*,*) 'Error setting an option'
  goto 9000

end program problemPCBeam
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!

subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt
  use dimopt
  implicit none
  integer N, NEW_X,I
  double precision F, X(N),sigmax(N),fmeantmp,fvartmp,fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision DAT(*)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision  rho, L, sigmay, pi, p, E, Fs 
  integer::myflag(10) 
        
  kprob=IDAT(1)
  probtype(1:N)=IDAT(3:N+2)

  do i=1,N
     sigmax(i)=DAT(i)
     Xsave(i)=X(i)
  end do

  !---- MEAN and VARIANCE OF worst OBJECTIVE FUNCTION

  
  if (IDAT(2).eq.1) then ! Deterministic with PC

     call CalcExact(X,N,fmeantmp,fmeanprimetmp,0,DAT(1001:1020))

     fvartmp=0.0d0
     fvarprimetmp=0.0d0

  else

     call PCestimate(N-2,N,x,sigmax,23,0,DAT(1001:1020),2,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
     if (fvartmp.lt.0.0) fvartmp=0.0
     
  end if
  
  do i=1,n
     data(i)=x(i)
  end do

  !---- COMBINED OBJECTIVE FUNCTION
  F=fmeantmp+fvartmp

  DAT(N+1)=fmeantmp
  DAT(N+2)=fvartmp

  do i=1,N
     DAT(N+2+i)= fmeanprimetmp(i)+fvarprimetmp(i)
  end do

 
  if (id_proc.eq.0) then
     print*,''
     write(*,'(4x,a,3F13.4)') '>>Objective:',fmeantmp,fvartmp,fmeantmp+fvartmp
     write(*,'(4x,a,F13.4)') '>>Coeff of variance :',sqrt(fvartmp)/fmeantmp
  end if

!  stop
  do i=1,n
     DAT(2*N+2+i)=Xsave(i)
     X(i)=Xsave(i)
  end do

  IERR = 0
  return

end subroutine EV_F

!
! =============================================================================
!
!                     Computation of constraints
!
! =============================================================================
!
subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt

  implicit none
  integer N, NEW_X, M
  double precision G(M), X(N), sigmax(N), cmean(M), cstd(M), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n),dc(M,N)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR, i, j, cnt
  double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision  rho, L, sigmay, pi, p, E, Fs 
  integer::myflag(10) 
      

  kprob=IDAT(1)
  probtype(1:N)=IDAT(3:N+2)
  
  do i=1,N
     sigmax(i)=DAT(i)
     Xsave(i)=X(i)
  end do

  dc(:,:)=0.0

  do i=1,M

     !---- MEAN OF INEQUALITY CONSTRAINT i


     if (IDAT(2).eq.1) then

        call CalcExact(X,N,fmeantmp,fmeanprimetmp,i,DAT(1001:1020))
        fvartmp=0.0
        fvarprimetmp(:)=0.0

     else

     call PCestimate(N-2,N,x,sigmax,23,i,DAT(1001:1020),2,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
     if (fvartmp.lt.0.0) fvartmp=0.0

     end if

!     if(id_proc.eq.0)   print*,"me,st:",fmeantmp,fvartmp

     cmean(i)=fmeantmp
     cstd(i)=sqrt(fvartmp)
     
!     if(id_proc.eq.0)   print*,"me,st:",cmean(i),cstd(i)
!     if(id_proc.eq.0)   print*,"mp,vp:",fmeanprimetmp(1:N),fvarprimetmp(1:N)


     do j=1,N
        dc(i,j)=fmeanprimetmp(j)
        if (fvartmp.ne.0.0) then
           dc(i,j)=dc(i,j)+kprob*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
        end if
     end do
  end do ! M Loop

  !---- COMBINED INEQUALITY CONSTRAINTS

  G(1:M)=cmean(1:M)+kprob*cstd(1:M)

  !Just printing

  if (id_proc.eq.0) then
     print*,''
     write(*,'(4x,a)') '>>Normalized Constraint Values:'
     do i=1,M
        write(*,'(3E13.5,f13.5)'),cmean(i),cstd(i),g(i),cmean(i)/cstd(i)
        DAT(1020+i)=(cmean(i)/cstd(i))
     end do
     print*,''
  end if

  !---- INEQUALITY CONSTRAINTS gradient

  do i=1,N
     DAT(3*N+2+i)=Xsave(i)
     X(i)=Xsave(i)
  end do

  cnt=0
  do i=1,M
     do j=1,N
        cnt=cnt+1
        DAT(4*N+2+cnt)=dc(i,j)
     end do
  end do

  IERR = 0
  return
end subroutine EV_G

!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt

  implicit none
  integer N, NEW_X,i
  double precision GRAD(N), X(N), sigmax(N), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision  rho, L, sigmay, pi, p, E, Fs 
  logical samex
  integer::myflag(10) 

  samex=.true.
  do i=1,N
     if (x(i).ne.DAT(2*N+2+i)) samex=.false. 
  end do
  
  if (samex) then

     !         if (id_proc.eq.0) print *,'Samex in obj',X

     !---- TOTAL GRADIENT OF OBJECTIVE FUNCTION

     do i=1,n
        GRAD(i)=DAT(N+2+i)
     end do

  else

     ! if (id_proc.eq.0) print *,'Not Samex in obj',X

     kprob=IDAT(1)
     probtype(1:N)=IDAT(3:N+2)

     do i=1,N
        sigmax(i)=DAT(i)
     end do

!!$      gtol=1e-4
!!$
!!$      low(1:N-1)=X(1:N-1)-sigmax(1:N-1)
!!$      up(1:N-1)=X(1:N-1)+sigmax(1:N-1)
!!$
!!$      call optimize(N-1,X,N,fmax,gradmax,low,up,gtol,.true.,.false.,1)


     !---- MEAN and VARIANCE OF worst OBJECTIVE FUNCTION



     if (IDAT(2).eq.1) then

        call CalcExact(X,N,fmeantmp,fmeanprimetmp,0,DAT(1001:1020))
        fvartmp=0.0
        fvarprimetmp(:)=0.0

     else

     call PCestimate(N-2,N,x,sigmax,23,0,DAT(1001:1020),2,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
     if (fvartmp.lt.0.0) fvartmp=0.0

     end if

     !---- OBJECTIVE FUNCTION gradient and x value

     do i=1,N
        grad(i)=fmeanprimetmp(i)+fvarprimetmp(i)
     end do
!!$         if (id_proc.eq.0) print *,'Obj Gradient',GRAD(1:3)
  end if

  IERR = 0
  return
end subroutine EV_GRAD_F

!
! =============================================================================
!
!                Computation of Jacobian of constraints
!
! =============================================================================
!
subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt
  
  implicit none
  integer TASK, N, NEW_X, M, NZ
  double precision X(N), A(NZ),dc(M,N), sigmax(N), fmeantmp, fvartmp
  integer ACON(NZ), AVAR(NZ), I, J, K, cnt, NMC
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision  rho, L, sigmay, pi, p, E, Fs
  integer IDAT(*)
  integer IERR, kprob
  logical samex
  integer::myflag(10) 

  if( TASK.eq.0 ) then 

     !
     !     structure of Jacobian:
     !
     
     
     ACON(1) = 1
     AVAR(1) = 1

     ACON(2) = 1
     AVAR(2) = 2

     ACON(3) = 1
     AVAR(3) = 3

     ACON(4) = 1
     AVAR(4) = 4


     ! con 2

     ACON(5) = 2
     AVAR(5) = 1

     ACON(6) = 2
     AVAR(6) = 2

     ACON(7) = 2
     AVAR(7) = 3

     ACON(8) = 2
     AVAR(8) = 4

     !con 3


     ACON(9)  = 3
     AVAR(9)  = 1

     ACON(10) = 3
     AVAR(10) = 2

     ACON(11) = 3
     AVAR(11) = 3

     ACON(12) = 3
     AVAR(12) = 4


  else

     samex=.true.
     do i=1,N
        if (x(i).ne.DAT(3*N+2+i)) samex=.false. 
     end do

     if (samex) then

        !            if (id_proc.eq.0) print *,'Samex in con'

        cnt=0
        do i=1,M
           do j=1,N
              cnt=cnt+1
              dc(i,j)=DAT(4*N+2+cnt)
           end do
        end do


     else

        !            if (id_proc.eq.0) print *,'Not Samex in con'

        !---- TOTAL GRADIENT OF CONSTRAINTS 

        kprob=IDAT(1)
        probtype(1:N)=IDAT(3:N+2)

        do i=1,N
           sigmax(i)=DAT(i)
        end do

        dc(:,:)=0.0

        do i=1,M

           !---- MEAN OF INEQUALITY CONSTRAINT i
           
           
           if (IDAT(2).eq.1) then
           
              call CalcExact(X,N,fmeantmp,fmeanprimetmp,i,DAT(1001:1020))
              fvartmp=0.0
              fvarprimetmp(:)=0.0

           else

              call PCestimate(N-2,N,x,sigmax,23,i,DAT(1001:1020),2,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
              if (fvartmp.lt.0.0) fvartmp=0.0

           end if

           do j=1,N
              dc(i,j)=fmeanprimetmp(j)
              if (fvartmp.ne.0.0) then
                 dc(i,j)=dc(i,j)+kprob*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
              endif
           end do

        end do

     end if
     
     ! Assemble

     A(1)=dc(1,1)
     A(2)=dc(1,2)
     A(3)=dc(1,3)
     A(4)=dc(1,4)


     A(5)=dc(2,1)
     A(6)=dc(2,2)
     A(7)=dc(2,3)
     A(8)=dc(2,4)

     A(9)=dc(3,1)
     A(10)=dc(3,2)
     A(11)=dc(3,3)
     A(12)=dc(3,4)

     !if (id_proc.eq.0) print *,'Cons Gradients',jac(1:6)
     
  end if


  IERR = 0
  return
end subroutine EV_JAC_G
!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
  implicit none
  integer TASK, N, NEW_X, M, NEW_LAM, NNZH,  i,j,ii
  double precision X(N), OBJFACT, LAM(M), HESS(NNZH), sigmax(N)
  integer IRNH(NNZH), ICNH(NNZH)
  double precision::fmeantmp,fvartmp
  double precision OBJHESS(NNZH),CONHESS(M,NNZH)
  double precision DAT(*)
  integer IDAT(*), kprob
  integer IERR
  integer::myflag(10) 
  
  if( TASK.eq.0 ) then

     !
     !     structure of sparse Hessian (lower triangle):
     !

     IRNH(1) = 1
     ICNH(1) = 1

     IRNH(2) = 2
     ICNH(2) = 2

     IRNH(3)=  3
     ICNH(3)=  3

     IRNH(4)=4
     ICNH(4)=4

     IRNH(5)= 2
     ICNH(5)=1


     IRNH(6)=3
     ICNH(6)=2

     IRNH(7)=4
     ICNH(7)=3

     IRNH(8)=3
     ICNH(8)=1

     IRNH(9)=4
     ICNH(9)=2

     IRNH(10) = 4
     ICNH(10) = 1


 
  else

     
     IERR = 1

  endif

  return
end subroutine EV_HESS



!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
  use dimopt
  use dimpce,only:probtype,id_proc,fcnt,fgcnt,fghcnt
  implicit none
  integer ALG_MODE, ITER_COUNT, LS_TRIAL
  double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
  double precision ALPHA_DU, ALPHA_PR
  double precision DAT(*),tol
  integer IDAT(*)
  integer ISTOP

  !
  !     You can put some output here
  !
  if (id_proc.eq.0) then

     if (ITER_COUNT .eq.0) then
        write(*,*) 
        write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
        write(86,*) 'iter    objective    betag1    betag2    betag3'
        write(37,*) 'iter    X1    X2        X3    X4'
                     
     end if

     write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU
     write(76,'(i5,5e15.7,3i8)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU,fcnt,fgcnt,fghcnt
     write(86,'(i5,4e15.7)') ITER_COUNT,OBJVAL,DAT(1020+1),DAT(1020+2),DAT(1020+3)
     write(37,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DATA(1),DATA(2),DATA(3),DATA(4)
 
  end if

  !
  !     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
  !     simple example.
  !
  
  if (ITER_COUNT .gt. 1 ) then
     
     open(unit=59,file='dnorm.inp',status='old')
     read(59,*) tol
     close(59)

     if (DNORM.le.tol) ISTOP = 1  

  end if



  
  return
end subroutine ITER_CB

!===================

subroutine epigrads(fct,fctindx,dim,ndimt,xtmp,xstdt,ftmp,dftmp)
  use omp_lib

  implicit none
  integer :: DIM,ndimt,fct,fctindx
  real*8,intent(in)  :: xtmp(ndimt),xstdt(ndimt)
  real*8,intent(out) :: dftmp(ndimt)
  real*8::ftmp
  real*8 :: gtol,low(ndimt-DIM),up(ndimt-DIM)


  gtol=1e-6

  low(1:ndimt-DIM)= xtmp(1:ndimt-DIM)  !+ xstdt(1:ndimt-DIM)
  up(1:ndimt-DIM) = xtmp(1:ndimt-DIM)  !+ xstdt(1:ndimt-DIM)

  call optimize(ndimt-DIM,xtmp,ndimt,ftmp,dftmp,low,up,gtol,.true.,.false.,fctindx)

!  call optimize(ndvar,D,ndvart,fobj,dfdD,low,up,gtol,maximize,outputscreen,fct)

  return
end subroutine epigrads

