  subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct)
    use dimpce,only:fctindx
    implicit none

    integer  :: ndvart,fct
    double precision :: X(ndvart),fobj,dfdD(ndvart),x3
    double precision ::  rho, L, sigmay, pi, p, E, Fs  
    
    fctindx=fct

    call my_calcf(x,ndvart,fobj)

    call my_calcdf(x,ndvart,dfdD)   

    return
  end subroutine CalcstuffBFGS
  
  
  !++++++++++++++++++++++++++++++++++

  subroutine CalcExact(X,ndvart,fobj,dfdD,fct,DATA)
    use dimpce,only:fctindx,DAT,mainprog
    implicit none

    real*8::DATA(20)
    integer  :: ndvart,fct
    double precision :: X(ndvart),fobj,dfdD(ndvart),x3
    double precision ::  rho, L, sigmay, pi, p, E, Fs  

    fctindx=fct
    mainprog=.false.
    DAT(1:20)=DATA(1:20)

    call my_calcf(x,ndvart,fobj)
    call my_calcdf(x,ndvart,dfdD) 

    return
  end subroutine CalcExact

  !+++++++++++++++++++++++++++++++++

  subroutine my_calcf(x,ndimt,f)

    use dimpce,only:fcnt,DAT,fctindx
    implicit none
    integer::ndimt
    real*8::x(ndimt),f
    real*8::sigma_allow,tau_allow,Fs,B,D,M,V


    fcnt=fcnt+1 !function counter

    ! reallocate supplied data

    sigma_allow=dat(1)
    tau_allow=dat(2)
    Fs=dat(3)

    B=x(1)
    D=x(2)
    M=x(3)
    V=x(4)

    if (fctindx.eq.0) then

       !---- OBJECTIVE FUNCTION
       f = B*D

    else if (fctindx.eq.1) then

       !---- INEQUALITY CONSTRAINTS
       !bending stress constraint

       f=(6.0*M*fs)/(b*(d**2)*sigma_allow)-1.0       

    else if (fctindx.eq.2) then

       ! Inequality constraint 2
       ! Shear Stress constraint

       f=(3.0*V*fs)/(2.0*b*d*tau_allow) -1.0

    else if (fctindx.eq.3) then

       f= d/(2.0*b) - 1.0

    else 

       print*, 'Wrong function index for this test case',fctindx
       stop
    end if


    return
  end subroutine my_calcf
  
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  subroutine my_calcdf(x,ndimt,df)

    use dimpce,only:fgcnt,DAT,fctindx
    implicit none
    integer::ndimt
    real*8::x(ndimt),df(ndimt)
    real*8::sigma_allow,tau_allow,Fs,B,D,M,V


    fgcnt=fgcnt+1 !function counter

    ! reallocate supplied data

    sigma_allow=dat(1)
    tau_allow=dat(2)
    Fs=dat(3)

    B=x(1)
    D=x(2)
    M=x(3)
    V=x(4)

    if (fctindx.eq.0) then

       !---- OBJECTIVE FUNCTION
       df(1)=d
       df(2)=b

    else if (fctindx.eq.1) then

       !---- INEQUALITY CONSTRAINT 1

       df(1)= -(6.0*fs*M)/((b*d)**2*sigma_allow)
       df(2)= -(12.0*M*fs)/(b*(d**3)*sigma_allow)


    else if (fctindx.eq.2) then

       !---- INEQUALITY CONSTRAINT 2

       df(1)=-(3.0*V*fs)/(2.0*(b**2)*d*tau_allow)
       df(2)=-(3.0*V*FS)/(2.0*b*(d**2)*tau_allow)

    else if (fctindx.eq.3) then

       df(1) = -1.0*d/(2.0*b**2)
       df(2) = 1.0/(2.0*b) 
    else

       print*, 'Wrong function index for this test case',fctindx
       stop

    end if


    return
  end subroutine my_calcdf
