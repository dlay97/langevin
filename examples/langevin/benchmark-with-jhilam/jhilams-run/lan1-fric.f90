!############ Langevin program for SF ##############
!##### 2D model in (Q_{20},Q_{30}) space ###########
!##### F90 implemented on Feb 7,21 #################

! Presumably pulls from section 2.3.2 of
! https://www.sciencedirect.com/science/article/abs/pii/0370157396000038

!===================================================================================================
!   Declare global variables
!===================================================================================================
module varbl
    implicit real*8 (a-h,o-z)
    implicit integer (i-n)

    double precision, dimension(:,:), allocatable :: &
        aa1,aa1c,aa1a,aa1cc,aa1aa,aa1ca,&
        aa2,aa2c,aa2a,aa2cc,aa2aa,aa2ca,&
        aa3,aa3c,aa3a,aa3cc,aa3aa,aa3ca

    double precision, dimension(:,:), allocatable :: &
        ee1,ee1c,ee1a,ee2,ee2c,ee2a,ee3,ee3c,ee3a,&
        gg1,gg1c,gg1a,gg2,gg2c,gg2a,gg3,gg3c,gg3a

    double precision, dimension(:,:), allocatable :: &
        ca,ff1,ff1c,ff1a,ff1cc,ff1aa,ff1ca,amas,charg

    integer, dimension(:,:), allocatable :: &
        isci

    double precision, dimension(:), allocatable :: &
        cc,aa

    double precision, dimension(:), allocatable :: &
        pro,cp,ap

    double precision, dimension(:), allocatable :: &
        cin,ain,pin

    integer icmax,iamax,npro
    double precision qci,qai,dcc,daa

end module varbl

!===================================================================================================
!   Main
!===================================================================================================
program main
    use varbl
    implicit real*8 (a-h,o-z)
    real ran
    external ran
    dimension x(10001),y(10001)
    dimension counta1(100),countz1(100),countzn(100,100)
    logical scin

!###########  input files ########################
    open(10,file='lan.in',status='old')
    open(20,file='prob.in',status='old')
    open(30,file='pot2030.in',status='old')
    open(40,file='mas2030.in',status='old')
    open(50,file='scission.in',status='old')
!########### output file ########################
    open(5,file='lan.out',status='unknown')

    num   =1000
    bgauss=0.25d0
    ! hbar  =6.4655d0
    aldf  =10.0d0
    ! hbar2 =hbar*hbar
    scin  =.false.

!########## reading from lan.in #################
    read(10,*)amass
    read(10,*)tmax,delt,e0
    read(10,*)nfric,diss11,diss22,diss12
    read(10,*)njob,npro
    read(10,*)qci,qcf,dcc
    read(10,*)qai,qaf,daa
    read(10,*)cstart,qsci
    read(10,*)kotl,alfa1,alfa2,alfstp
    read(10,*)ia1,ia2
    read(10,*)iz1,iz2
    read(10,*)in1,in2
!###### writing Langevin parameters ############
    write(5,*)'CN mass =', amass
    write(5,*)'max Langevin time =',tmax
    write(5,*)'time step size =',delt
    write(5,*)'value of E0 (MeV) =',e0
    write(5,*)'range of Q20 for input data from ',qci,' to ',qcf,' in steps of ',dcc
    write(5,*)'range of Q30 for input data from ',qai,' to ',qaf,' in steps of ',daa
    Write(5,*)'range of coordinates for calculations'
    write(5,*)'starting ref. for Q20 =',cstart
    if(kotl.eq.1)then
        write(5,*)'Q30 values on OTL from',alfa11,'to',alfa21,'in steps of',alfstp1
    else
        write(5,*)'Q20-Q30 values on OTL are same as WKB-prob input file'
    endif
    write(5,*)'Q_N value for scission condition = ',qsci

    icmax=(qcf-qci)/dcc+1.0d0
    iamax=(qaf-qai)/daa+1.0d0
!######## allocate all variables ###############
    allocate(  aa1(icmax,iamax))
    allocate( aa1c(icmax,iamax))
    allocate( aa1a(icmax,iamax))
    allocate(aa1cc(icmax,iamax))
    allocate(aa1aa(icmax,iamax))
    allocate(aa1ca(icmax,iamax))
    allocate(  aa2(icmax,iamax))
    allocate( aa2c(icmax,iamax))
    allocate( aa2a(icmax,iamax))
    allocate(aa2cc(icmax,iamax))
    allocate(aa2aa(icmax,iamax))
    allocate(aa2ca(icmax,iamax))
    allocate(  aa3(icmax,iamax))
    allocate( aa3c(icmax,iamax))
    allocate( aa3a(icmax,iamax))
    allocate(aa3cc(icmax,iamax))
    allocate(aa3aa(icmax,iamax))
    allocate(aa3ca(icmax,iamax))

    allocate(  ee1(icmax,iamax))
    allocate( ee1c(icmax,iamax))
    allocate( ee1a(icmax,iamax))
    allocate(  ee2(icmax,iamax))
    allocate( ee2c(icmax,iamax))
    allocate( ee2a(icmax,iamax))
    allocate(  ee3(icmax,iamax))
    allocate( ee3c(icmax,iamax))
    allocate( ee3a(icmax,iamax))
    allocate(  gg1(icmax,iamax))
    allocate( gg1c(icmax,iamax))
    allocate( gg1a(icmax,iamax))
    allocate(  gg2(icmax,iamax))
    allocate( gg2c(icmax,iamax))
    allocate( gg2a(icmax,iamax))
    allocate(  gg3(icmax,iamax))
    allocate( gg3c(icmax,iamax))
    allocate( gg3a(icmax,iamax))

    allocate(  ff1(icmax,iamax))
    allocate( ff1c(icmax,iamax))
    allocate( ff1a(icmax,iamax))
    allocate(ff1cc(icmax,iamax))
    allocate(ff1aa(icmax,iamax))
    allocate(ff1ca(icmax,iamax))
    allocate( isci(icmax,iamax))
    allocate( amas(icmax,iamax))
    allocate(charg(icmax,iamax))

    allocate(   cc(icmax))
    allocate(   ca(icmax,iamax+1))
    allocate(   aa(iamax))
    allocate(   pro(npro))
    allocate(    cp(npro))
    allocate(    ap(npro))

!######## reading WKB probabilities & OTL ######
    do i=1,npro
        read(20,*)cp(i),ap(i),pro(i)
    enddo

!######## reads potential energy ###############
    do ic=1,icmax
        do ia=1,iamax
            read(30,*)cc(ic),aa(ia),ff1(ic,ia)
            ca(ic,1)=cc(ic)
            ca(ic,ia+1)=aa(ia)
        enddo
    enddo
!######## reads inertia tensor #################
    do ic=1,icmax
        do ia=1,iamax
            read(40,*)cdum,adum,am1,am2,am3
!   am1=1.0d0
!   am2=0.0d0
!   am3=1.0d0
!######## inverting mass tensor ################
            det=am1*am3-am2*am2
            if(det.lt.0.000001d0)det=1.0d0
            aa1(ic,ia)= am3/det
            aa2(ic,ia)=-am2/det
            aa3(ic,ia)= am1/det
        enddo
    enddo
    
!######## reads scission line & fragments ######
    do ic=1,icmax
        do ia=1,iamax
            read(50,*)cdum,adum,qn,charg(ic,ia),amas(ic,ia)
            if(qn.gt.qsci)then
                isci(ic,ia)=0
            else
                isci(ic,ia)=1
            endif
        enddo
    enddo
    close(10)
    close(20)
    close(30)
    close(40)
    close(50)
!###### all reading done #######################

!######## printing the scission line ###########
! write(5,*)'coordiantes of scission line'
! do i=1,icmax
!  do j=1,iamax
!   if((isci(i,j).eq.1.and.isci(i,j+1).eq.0).or. &
!      (isci(i,j).eq.0.and.isci(i,j+1).eq.1))then
!    write(5,32)cc(i),aa(j)
!   endif
!  enddo
! enddo
!32 format(2f10.2)

!## calculate derivatives of inverse mass tensor and PES
    call prep_tensor_derivative_dx(aa1,aa1c)
    call prep_tensor_derivative_dy(aa1,aa1a)
    call prep_tensor_second_derivatives (aa1c,aa1a,aa1cc,aa1aa,aa1ca)
    call prep_tensor_derivative_dx(aa2,aa2c)
    call prep_tensor_derivative_dy(aa2,aa2a)
    call prep_tensor_second_derivatives (aa2c,aa2a,aa2cc,aa2aa,aa2ca)
    call prep_tensor_derivative_dx(aa3,aa3c)
    call prep_tensor_derivative_dy(aa3,aa3a)
    call prep_tensor_second_derivatives (aa3c,aa3a,aa3cc,aa3aa,aa3ca)
    call prep_tensor_derivative_dx(ff1,ff1c)
    call prep_tensor_derivative_dy(ff1,ff1a)
    call prep_tensor_second_derivatives (ff1c,ff1a,ff1cc,ff1aa,ff1ca)
!## define friction tensor #####################
    do ic=1,icmax
        do ia=1,iamax
            ee1(ic,ia)=diss11
            ee3(ic,ia)=diss22
            ee2(ic,ia)=diss12
        enddo
    enddo
!## derivatives of friction tensor #############
    call prep_tensor_derivative_dx(ee1,ee1c)
    call prep_tensor_derivative_dy(ee1,ee1a)
    call prep_tensor_derivative_dx(ee2,ee2c)
    call prep_tensor_derivative_dy(ee2,ee2a)
    call prep_tensor_derivative_dx(ee3,ee3c)
    call prep_tensor_derivative_dy(ee3,ee3a)
    ! write(*,*) ee1c
    ! stop
!## random force strength tensor ###############
    do ic=1,icmax
        do ia=1,iamax
            et1=ee1(ic,ia)
            et2=ee2(ic,ia)
            et3=ee3(ic,ia)
            call get_matrix_square_root(ci,ai,et1,et2,et3,r1,r2,r3)
            gg1(ic,ia)=r1
            gg2(ic,ia)=r2
            gg3(ic,ia)=r3
        enddo
    enddo
!  stop
!## derivatives of random force strength #######
    call prep_tensor_derivative_dx(gg1,gg1c)
    call prep_tensor_derivative_dy(gg1,gg1a)
    call prep_tensor_derivative_dx(gg2,gg2c)
    call prep_tensor_derivative_dy(gg2,gg2a)
    call prep_tensor_derivative_dx(gg3,gg3c)
    call prep_tensor_derivative_dy(gg3,gg3a)
!## level density parameter and initial e* #####
    alvl=amass/aldf
    exint=0.0d0
    temp=0.0d0

!## table for gaussian random no. ##############
    call table1(bgauss,x,y,dx,num)

!## setting the outer turning points ###########
    call set_langevin_starting_points(kotl,alfa1,alfa2,alfstp,cstart,ki)
    ! write(*,*) kotl
    ! write(*,*) ki
    
    ! kotl, alfa1, alfa2, alfstep, cstart
    ! do i1=1,

!---initial coordinates and momenta-------
    counta1=0.0d0
    countz1=0.0d0
    countzn=0.0d0
    nna1=(ia2-ia1)+1
    nnz1=(iz2-iz1)+1
    nnzn=(in2-in1)+1
    akk=0.0d0
    kke=0
    delt2=dsqrt(delt)
    delt3=delt*delt2
    deltt=0.5*delt*delt

!############ loop for different points on the OTL #############
    do it=1,ki
!############ loop for no. of events for a particulat OTL ######
        do irun=1,njob
            call CPU_TIME(startTime)
!############ initialization for an event ######################
            pp0=pin(it)
            qc0=cin(it)
            qa0=ain(it)
            qc1=qc0
            qa1=qa0
            pc1=0.0d0
            pa1=0.0d0
            l=0
            trun=0.0d0
            nofric=1
!############ time propagation starts ###########################
            do
!----------------------------------------------------------------
                call interpolate_pes_and_derivatives(qc1,qa1,f1,f1c,f1a,f1cc,f1aa,f1ca)
                ! write(*,*) qc1, qa1
                ! write(*,*) f1,f1c,f1a,f1cc,f1aa,f1ca
                ! stop
!----------------------------------------------------------------
                call interpolate_tensors_and_derivatives(qc1,qa1,a1,a2,a3,e1,e2,e3,g1,g2,g3, &
                    a1c,a2c,a3c,a1a,a2a,a3a,a1cc,a2cc,a3cc,a1aa,a2aa,a3aa,&
                    a1ca,a2ca,a3ca,e1c,e2c,e3c,e1a,e2a,e3a,g1c,g2c,g3c,g1a,g2a,g3a)
                
!----------------------------------------------------------------
                exint=e0-0.5d0*pc1*pc1*a1-0.5d0*pa1*pa1*a3-pc1*pa1*a2-f1
                
!    write(*,*)qc1,qa1,e0,exint,f1,trun

!---selecting the dynamics with or without friction--------------
                if(nfric.eq.0)then
                    nofric=1
                    goto 654
                endif
                if(exint.le.0.0d0)then
                    nofric=1
                else
                    nofric=0
                endif
                if(nofric.eq.1)goto 654
!----------------------------------------------------------------
                !=====================================================================================
                !
                ! Basically everything up to the scission check follows from
                ! https://www.sciencedirect.com/science/article/abs/pii/0370157396000038 Appendix A.2
                !
                !=====================================================================================

                temp=1.!dsqrt(exint/alvl)
                call gran(idum,rand1,x,y,dx,num)
                call gran(idum,rand2,x,y,dx,num)
                call gran(idum,rand3,x,y,dx,num)
                call gran(idum,rand4,x,y,dx,num)
                ! rand1 = 1.
                ! rand2 = 1.
                ! rand3 = 1.
                ! rand4 = 1.
                ! write(*,*) rand1, rand2, rand3, rand4
                ! nofric = 0
                ! write(*,*) 'nofric',nofric
                t1=dsqrt(temp)
                g1 =g1 *t1
                g2 =g2 *t1
                g3 =g3 *t1
                g1c=g1c*t1
                g2c=g2c*t1
                g3c=g3c*t1
                g1a=g1a*t1
                g2a=g2a*t1
                g3a=g3a*t1
!---------------------------------------------------------------
                ex1=(e1*a1+e2*a2)
                ex2=(e1*a2+e2*a3)
                ex3=(e3*a3+e2*a2)
                ex4=(e3*a2+e2*a1)
                r1c =delt2*rand1
                r1a =delt2*rand2
                r2c =delt3*(0.5d0*rand1+0.288675d0*rand3)
                r2a =delt3*(0.5d0*rand2+0.288675d0*rand4)
                rr2c=delt3*(0.5d0*rand1-0.288675d0*rand3)
                rr2a=delt3*(0.5d0*rand2-0.288675d0*rand4)
!--------------------------------------------------------------
654             continue
                vc1  =pc1*a1+pa1*a2
                va1  =pa1*a3+pc1*a2
                vc1c =pc1*a1c+pa1*a2c
                vc1a =pc1*a1a+pa1*a2a
                vc1pc=a1
                vc1pa=a2
                va1c =pa1*a3c+pc1*a2c
                va1a =pa1*a3a+pc1*a2a
                va1pc=a2
                va1pa=a3
                ! write(*,*) f1,f1c,f1a,f1cc,f1aa,f1ca
                ! write(*,*) 'mass'
                ! write(*,*) a1,a2,a3
                ! write(*,*) a1c,a1a, a2c,a2a, a3c, a3a
                ! write(*,*) a1cc, a1ca, a1aa, a2cc, a2ca, a2aa, a3cc, a3ca, a3aa
                ! write(*,*) 'v',vc1,va1
                ! write(*,*) 'fric',e1,e2,e3
!--------------------------------------------------------------
                if(nofric.eq.1)then
                    h1=-f1c-0.5d0*a1c*pc1*pc1-0.5d0*a3c*pa1*pa1-a2c*pc1*pa1
                    h1c=-f1cc-0.5d0*a1cc*pc1*pc1-0.5d0*a3cc*pa1*pa1-a2cc*pc1*pa1
                    h1a=-f1ca-0.5d0*a1ca*pc1*pc1-0.5d0*a3ca*pa1*pa1-a2ca*pc1*pa1
                    h1pc=-a1c*pc1-a2c*pa1
                    h1pa=-a3c*pa1-a2c*pc1

                    h2=-f1a-0.5d0*a1a*pc1*pc1-0.5d0*a3a*pa1*pa1-a2a*pc1*pa1
                    h2c=-f1ca-0.5d0*a1ca*pc1*pc1-0.5d0*a3ca*pa1*pa1-a2ca*pc1*pa1
                    h2a=-f1aa-0.5d0*a1aa*pc1*pc1-0.5d0*a3aa*pa1*pa1-a2aa*pc1*pa1
                    h2pc=-a1a*pc1-a2a*pa1
                    h2pa=-a3a*pa1-a2a*pc1

                    pc2=pc1+h1*delt+deltt*(vc1*h1c+va1*h1a+h1*h1pc+h2*h1pa)
                    pa2=pa1+h2*delt+deltt*(vc1*h2c+va1*h2a+h1*h2pc+h2*h2pa)

                    qc2=qc1+vc1*delt+deltt*(vc1*vc1c+va1*vc1a+h1*vc1pc+h2*vc1pa)
                    qa2=qa1+va1*delt+deltt*(vc1*va1c+va1*va1a+h1*va1pc+h2*va1pa)
!--------------------------------------------------------------
                else
!--------------------------------------------------------------
                    h1=-f1c-0.5d0*a1c*pc1*pc1-0.5d0*a3c*pa1*pa1-a2c*pc1*pa1 &
                        -ex1*pc1-ex2*pa1 !issue here: should contract a1c with coordinates, not momentum

                    h1c=-f1cc-0.5d0*a1cc*pc1*pc1-0.5d0*a3cc*pa1*pa1-a2cc*pc1*pa1 &
                        -(e1c*a1+e1*a1c+e2c*a2+e2*a2c)*pc1 &
                        -(e1c*a2+e1*a2c+e2c*a3+e2*a3c)*pa1

                    h1a=-f1ca-0.5d0*a1ca*pc1*pc1-0.5d0*a3ca*pa1*pa1-a2ca*pc1*pa1 &
                        -(e1a*a1+e1*a1a+e2a*a2+e2*a2a)*pc1 &
                        -(e1a*a2+e1*a2a+e2a*a3+e2*a3a)*pa1

                    h1pc=-a1c*pc1-a2c*pa1-e1*a1-e2*a2
                    h1pa=-a3c*pa1-a2c*pc1-e1*a2-e2*a3

                    h2=-f1a-0.5d0*a1a*pc1*pc1-0.5d0*a3a*pa1*pa1-a2a*pc1*pa1 &
                        -ex3*pa1-ex4*pc1

                    h2c=-f1ca-0.5d0*a1ca*pc1*pc1-0.5d0*a3ca*pa1*pa1-a2ca*pc1*pa1 &
                        -(e3c*a3+e3*a3c+e2c*a2+e2*a2c)*pa1 &
                        -(e3c*a2+e3*a2c+e2c*a1+e2*a1c)*pc1

                    h2a=-f1aa-0.5d0*a1aa*pc1*pc1-0.5d0*a3aa*pa1*pa1-a2aa*pc1*pa1 &
                        -(e3a*a3+e3*a3a+e2a*a2+e2*a2a)*pa1 &
                        -(e3a*a2+e3*a2a+e2a*a1+e2*a1a)*pc1

                    h2pc=-a1a*pc1-a2a*pa1-e3*a2-e2*a1
                    h2pa=-a3a*pa1-a2a*pc1-e3*a3-e2*a2

                    ! write(*,*) '-f1c-0.5d0*a1c*pc1*pc1-0.5d0*a3c*pa1*pa1-a2c*pc1*pa1',-f1c-0.5d0*a1c*pc1*pc1&
                    !         -0.5d0*a3c*pa1*pa1-a2c*pc1*pa1
                    ! write(*,*) 'h',h1,h2
                    ! write(*,*) 'dhdq',h1c,h1a,h2c,h2a

                    pc2=pc1+h1*delt+deltt*(vc1*h1c+va1*h1a+h1*h1pc+h2*h1pa) &
                        +g1*r1c+g2*r1a+r2c*(g1*h1pc+g2*h1pa)+r2a*(g2*h1pc+g3*h1pa) &
                        +rr2c*(vc1*g1c+va1*g1a)+rr2a*(vc1*g2c+va1*g2a)

                    pa2=pa1+h2*delt+deltt*(vc1*h2c+va1*h2a+h1*h2pc+h2*h2pa) &
                        +g2*r1c+g3*r1a+r2c*(g1*h2pc+g2*h2pa)+r2a*(g2*h2pc+g3*h2pa) &
                        +rr2c*(vc1*g2c+va1*g2a)+rr2a*(vc1*g3c+va1*g3a)

                    qc2=qc1+vc1*delt+deltt*(vc1*vc1c+va1*vc1a+h1*vc1pc+h2*vc1pa) &
                        +r2c*(g1*vc1pc+g2*vc1pa)+r2a*(g2*vc1pc+g3*vc1pa)

                    qa2=qa1+va1*delt+deltt*(vc1*va1c+va1*va1a+h1*va1pc+h2*va1pa) &
                        +r2c*(g1*va1pc+g2*va1pa)+r2a*(g2*va1pc+g3*va1pa)
!-----------------------------------------------------------------
                endif
                ! write(*,*) 'pc2,pa2',pc2,pa2
                ! write(*,*) 'qc2,qa2',qc2,qa2
                ! stop
                l=l+1
                if(abs(qc1-qc2).gt.5.or.abs(qa1-qa2).gt.2.or.qc2.lt.cstart.or. &
                    qa2.gt.aa(iamax-2).or.qa2.lt.aa(1))exit
                pc1=pc2
                pa1=pa2
                qc1=qc2
                qa1=qa2
                ! write(*,*) 'l',l,'qc',qc1,'qa',qa1
!################### scission check ##############################
                ics=(qc1-ca(1,1))/dcc+1.0000000001d0
                ias=(qa1-ca(1,2))/daa+1.0000000001d0
                if(isci(ics,ias).eq.1.or.isci(ics+1,ias).eq.1.or. &
                    isci(ics,ias+1).eq.1.or.isci(ics+1,ias+1).eq.1)then
                    scin=.true.
                    exit
                endif
!################## max time check ###############################
                trun=trun+delt
                if(trun.ge.tmax)then
                    write(*,*)'trun gt tmax in ',irun,' of ',it,qc1,qa1
                    exit
                endif
            enddo
!########## frag distribution for scission #######################
            if(scin)then
                call inttwo(qc1,qa1,amas,fraga)
                call inttwo(qc1,qa1,charg,fragz)
                fragn=fraga-fragz
                akk=akk+pp0
                kke=kke+1
                do ii=1,nna1
                    ia=ia1+ii-1
                    if(nint(fraga).eq.ia)counta1(ii)=counta1(ii)+pp0
                enddo
                do ii=1,nnz1
                    iz=iz1+ii-1
                    if(nint(fragz).eq.iz) countz1(ii)=countz1(ii)+pp0
                enddo
                do ii=1,nnz1
                    iz=iz1+ii-1
                    do jj=1,nnzn
                        in=in1+jj-1
                        if(nint(fragz).eq.iz.and.nint(fragn).eq.in)then
                            countzn(ii,jj)=countzn(ii,jj)+pp0
                        endif
                    enddo
                enddo
            endif
            call CPU_TIME(endTime)
            write(*,*)'=============================================================='
            ! write(*,*) 'Endpoint number: ', it, 'Run number: ',irun,'Iteration time: ',endTime-startTime
!-------------------------------------------------------------------
        enddo
    enddo
    write(5,*) 'asdf'
if(kke.eq.0)then
 write(5,*)'no fission event'
else
  write(5,*)'****  mass distribution  ****'
  do ii=1,nna1
   prob1=counta1(ii)/akk
   write(5,23)ia1+ii-1, prob1
  enddo

  write(5,*)'****  charge distribution ****'
  do ii=1,nnz1
   probz1=countz1(ii)/akk
   write(5,23)iz1+ii-1,probz1
  enddo

  write(5,*)'Z-N correlation*******'
  do ii=1,nnz1
   do jj=1,nnzn
    probzn=countzn(ii,jj)/akk
    write(5,13)iz1+ii-1,in1+jj-1,probzn
   enddo
  enddo
 endif
 13   format(1x,2i10,f12.3)
 23   format(1x, i10,f12.3)
    stop
end

!===================================================================================================
!    Evaluating the potential energy and its first and second derivatives at point (arg1, arg2)
!===================================================================================================
subroutine interpolate_pes_and_derivatives(arg1,arg2,f1,f1c,f1a,f1cc,f1aa,f1ca)
    use varbl
    implicit real*8 (a-h,o-z)
    c=arg1
    alfa=arg2
    ic=(c-ca(1,1))/dcc+1.00000001d0
    ia=(alfa-ca(1,2))/daa+1.00000001d0
! write(*,*) 'arg1', arg1, 'arg2', arg2
!If on a boundary, return boundary value
    if(ic.gt.icmax-1.or.ia.gt.iamax-1)then
        if(ic.gt.icmax-1) then
            f1=ff1(icmax,ia)
            f1c=ff1c(icmax,ia)
            f1a=ff1a(icmax,ia)
            f1cc=ff1cc(icmax,ia)
            f1aa=ff1aa(icmax,ia)
            f1ca=ff1ca(icmax,ia)
        endif
        if(ia.gt.iamax-1) then
            f1=ff1(ic,iamax)
            f1c=ff1c(ic,iamax)
            f1a=ff1a(ic,iamax)
            f1cc=ff1cc(ic,iamax)
            f1aa=ff1aa(ic,iamax)
            f1ca=ff1ca(ic,iamax)
        endif
!csy  write(*,*)'error in interpolate_pes_and_derivatives  ','c=',c,'alfa=',alfa
        return
    endif

!Otherwise, call the interpolator
    p=(c-dcc*dfloat(ic-1)-qci)/dcc
    q=(alfa-daa*dfloat(ia-1)-qai)/daa
    
!-----calculates function values-------------------------
    call bilinear_interp(ic,ia,p,q,ff1,f1)
!----- calculates first derivatives----------------------
    call bilinear_interp(ic,ia,p,q,ff1c,f1c)
    call bilinear_interp(ic,ia,p,q,ff1a,f1a)
!----- calculates second derivatives---------------------
    call bilinear_interp(ic,ia,p,q,ff1cc,f1cc)
    call bilinear_interp(ic,ia,p,q,ff1aa,f1aa)
    call bilinear_interp(ic,ia,p,q,ff1ca,f1ca)
    return
end subroutine

!===================================================================================================
!    Evaluates the tensors aa, ee, and gg and their first and second derivatives at the point
!    (arg1, arg2)
!===================================================================================================
subroutine interpolate_tensors_and_derivatives(arg1,arg2,a1,a2,a3,e1,e2,e3,g1,g2,g3,a1c,a2c, &
    a3c,a1a,a2a,a3a,a1cc,a2cc,a3cc,a1aa,a2aa,a3aa,a1ca,a2ca,a3ca,e1c, &
    e2c,e3c,e1a,e2a,e3a,g1c,g2c,g3c,g1a,g2a,g3a)
    use varbl
    implicit real*8 (a-h,o-z)
    c=arg1
    alfa=arg2
    ic=(c-ca(1,1))/dcc+1.00000001d0
    ia=(alfa-ca(1,2))/daa+1.00000001d0
    if(ic.gt.icmax-1.or.ia.gt.iamax-1)then
        if(ic.gt.icmax-1) then
            a1=aa1(icmax,ia)
            a2=aa2(icmax,ia)
            a3=aa3(icmax,ia)
            e1=ee1(icmax,ia)
            e2=ee2(icmax,ia)
            e3=ee3(icmax,ia)
            g1=gg1(icmax,ia)
            g2=gg2(icmax,ia)
            g3=gg3(icmax,ia)
            a1c=aa1c(icmax,ia)
            a2c=aa2c(icmax,ia)
            a3c=aa3c(icmax,ia)
            a1a=aa1a(icmax,ia)
            a2a=aa2a(icmax,ia)
            a3a=aa3a(icmax,ia)
            e1c=ee1c(icmax,ia)
            e2c=ee2c(icmax,ia)
            e3c=ee3c(icmax,ia)
            e1a=ee1a(icmax,ia)
            e2a=ee2a(icmax,ia)
            e3a=ee3a(icmax,ia)
            g1c=gg1c(icmax,ia)
            g2c=gg2c(icmax,ia)
            g3c=gg3c(icmax,ia)
            g1a=gg1a(icmax,ia)
            g2a=gg2a(icmax,ia)
            g3a=gg3a(icmax,ia)
            a1cc=aa1cc(icmax,ia)
            a2cc=aa2cc(icmax,ia)
            a3cc=aa3cc(icmax,ia)
            a1aa=aa1aa(icmax,ia)
            a2aa=aa2aa(icmax,ia)
            a3aa=aa3aa(icmax,ia)
            a1ca=aa1ca(icmax,ia)
            a2ca=aa2ca(icmax,ia)
            a3ca=aa3ca(icmax,ia)
        endif
        if(ia.gt.iamax-1) then
            a1=aa1(ic,iamax)
            a2=aa2(ic,iamax)
            a3=aa3(ic,iamax)
            e1=ee1(ic,iamax)
            e2=ee2(ic,iamax)
            e3=ee3(ic,iamax)
            g1=gg1(ic,iamax)
            g2=gg2(ic,iamax)
            g3=gg3(ic,iamax)
            a1c=aa1c(ic,iamax)
            a2c=aa2c(ic,iamax)
            a3c=aa3c(ic,iamax)
            a1a=aa1a(ic,iamax)
            a2a=aa2a(ic,iamax)
            a3a=aa3a(ic,iamax)
            e1c=ee1c(ic,iamax)
            e2c=ee2c(ic,iamax)
            e3c=ee3c(ic,iamax)
            e1a=ee1a(ic,iamax)
            e2a=ee2a(ic,iamax)
            e3a=ee3a(ic,iamax)
            g1c=gg1c(ic,iamax)
            g2c=gg2c(ic,iamax)
            g3c=gg3c(ic,iamax)
            g1a=gg1a(ic,iamax)
            g2a=gg2a(ic,iamax)
            g3a=gg3a(ic,iamax)
            a1cc=aa1cc(ic,iamax)
            a2cc=aa2cc(ic,iamax)
            a3cc=aa3cc(ic,iamax)
            a1aa=aa1aa(ic,iamax)
            a2aa=aa2aa(ic,iamax)
            a3aa=aa3aa(ic,iamax)
            a1ca=aa1ca(ic,iamax)
            a2ca=aa2ca(ic,iamax)
            a3ca=aa3ca(ic,iamax)
        endif
!csy  write(*,*)'error in interpolate_tensors_and_derivatives  ','c=',c,'alfa=',alfa
        return
    endif
    p=(c-dcc*dfloat(ic-1)-qci)/dcc
    q=(alfa-daa*dfloat(ia-1)-qai)/daa
!-----calculates function values-------------------------
    call bilinear_interp(ic,ia,p,q,aa1,a1)
    call bilinear_interp(ic,ia,p,q,aa2,a2)
    call bilinear_interp(ic,ia,p,q,aa3,a3)
    call bilinear_interp(ic,ia,p,q,ee1,e1)
    call bilinear_interp(ic,ia,p,q,ee2,e2)
    call bilinear_interp(ic,ia,p,q,ee3,e3)
    call bilinear_interp(ic,ia,p,q,gg1,g1)
    call bilinear_interp(ic,ia,p,q,gg2,g2)
    call bilinear_interp(ic,ia,p,q,gg3,g3)
!----- calculates first derivatives------------------------
    call bilinear_interp(ic,ia,p,q,aa1c,a1c)
    call bilinear_interp(ic,ia,p,q,aa2c,a2c)
    call bilinear_interp(ic,ia,p,q,aa3c,a3c)
    call bilinear_interp(ic,ia,p,q,aa1a,a1a)
    call bilinear_interp(ic,ia,p,q,aa2a,a2a)
    call bilinear_interp(ic,ia,p,q,aa3a,a3a)
    call bilinear_interp(ic,ia,p,q,ee1c,e1c)
    call bilinear_interp(ic,ia,p,q,ee2c,e2c)
    call bilinear_interp(ic,ia,p,q,ee3c,e3c)
    call bilinear_interp(ic,ia,p,q,ee1a,e1a)
    call bilinear_interp(ic,ia,p,q,ee2a,e2a)
    call bilinear_interp(ic,ia,p,q,ee3a,e3a)
    call bilinear_interp(ic,ia,p,q,gg1c,g1c)
    call bilinear_interp(ic,ia,p,q,gg2c,g2c)
    call bilinear_interp(ic,ia,p,q,gg3c,g3c)
    call bilinear_interp(ic,ia,p,q,gg1a,g1a)
    call bilinear_interp(ic,ia,p,q,gg2a,g2a)
    call bilinear_interp(ic,ia,p,q,gg3a,g3a)
!----- calculates second derivatives------------------------
    call bilinear_interp(ic,ia,p,q,aa1cc,a1cc)
    call bilinear_interp(ic,ia,p,q,aa2cc,a2cc)
    call bilinear_interp(ic,ia,p,q,aa3cc,a3cc)
    call bilinear_interp(ic,ia,p,q,aa1aa,a1aa)
    call bilinear_interp(ic,ia,p,q,aa2aa,a2aa)
    call bilinear_interp(ic,ia,p,q,aa3aa,a3aa)
    call bilinear_interp(ic,ia,p,q,aa1ca,a1ca)
    call bilinear_interp(ic,ia,p,q,aa2ca,a2ca)
    call bilinear_interp(ic,ia,p,q,aa3ca,a3ca)
    return
end subroutine

!===================================================================================================
!    I think this is bilinear interpolation: https://en.wikipedia.org/wiki/Bilinear_interpolation
!    but for a weird data structure. There's more down in inttwo, where it handles cases where
!    we're off the grid
!===================================================================================================
subroutine bilinear_interp(ic,ia,p,q,quan,valu)
    use varbl
    implicit real*8 (a-h,o-z)
    dimension quan(icmax,iamax)

    f00=quan(ic,ia)
    f10=quan(ic+1,ia)
    f01=quan(ic,ia+1)
    f11=quan(ic+1,ia+1)

    valu=f00+p*(f10-f00)+q*(f01-f00)+p*q*(f00+f11-f01-f10)

    return
end subroutine

!--------------------------------------------------------
subroutine prep_tensor_derivative_dx(ff,dffdc)
    use varbl
!-- calculates derivative of ff with respect to c, extrapolated values used
!-- on boundary points
    implicit real*8 (a-h,o-z)
    dimension ff(icmax,iamax),dffdc(icmax,iamax)
    do i=2,icmax-1
        do j=1,iamax
            dffdc(i,j)=(ff(i+1,j)-ff(i-1,j))/(2.0d0*dcc)
        enddo
    enddo
    do j=1,iamax
        dffdc(1,j)=2.0d0*dffdc(2,j)-dffdc(3,j)
    enddo
    do j=1,iamax
        dffdc(icmax,j)=2.0d0*dffdc(icmax-1,j)-dffdc(icmax-2,j)
    enddo
    return
end subroutine

!--------------------------------------------------------
subroutine prep_tensor_derivative_dy(ff,dffda)
    use varbl
!-- calculates derivative of ff with respect to h, extrapolated values used
!-- on boundary points
    implicit real*8 (a-h,o-z)
    dimension ff(icmax,iamax),dffda(icmax,iamax)
    do i=1,icmax
        do j=2,iamax-1
            dffda(i,j)=(ff(i,j+1)-ff(i,j-1))/(2.0d0*daa)
        enddo
        dffda(i,1)=2.0d0*dffda(i,2)-dffda(i,3)
        dffda(i,iamax)=2.0d0*dffda(i,iamax-1)-dffda(i,iamax-2)
    enddo
    return
end subroutine

!--------------------------------------------------------
subroutine prep_tensor_second_derivatives(dffdc,dffda,d2fcc,d2faa,d2fca)
    use varbl
!-- must call prep_tensor_derivative_dx and prep_tensor_derivative_dy before calling this routine
!-- calculates second derivative of ff with respect to c,alfa,
!-----and (c,alfa),
!-- uses first derivative  values on the grid
    implicit real*8 (a-h,o-z)
    dimension dffdc(icmax,iamax),dffda(icmax,iamax)
    dimension d2fcc(icmax,iamax),d2faa(icmax,iamax),d2fca(icmax,iamax)
    call prep_tensor_derivative_dx(dffdc,d2fcc)
    call prep_tensor_derivative_dy(dffda,d2faa)
    call prep_tensor_derivative_dy(dffdc,d2fca)
    return
end subroutine

!===================================================================================================
!    Computes the square root of a matrix $f$, i.e. the symmetric matrix $r$ such that
!    $r_{ij}r_{jk}=f_{ik}$. A python equivalent is
!    https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.sqrtm.html
!===================================================================================================
subroutine get_matrix_square_root(c,alfa,fcc,fca,faa,rcc,rca,raa)
!-- calculates strength tensor of random force(rcc,rch,rhh), august 14, 2003
    implicit real*8(a-h,o-z)
    nz=100000
    eps=1.0d-10

!--following reads the friction tensor--------------------
    det=fcc*faa-fca**2
    if(det.le.0.0d0 ) then
        write(*,*)'friction in error for c,alfa=',c,alfa
    endif
    if(dabs(fcc).lt.eps.or.dabs(faa).lt.eps) then
        rcc=0.0d0
        rca=0.0d0
        raa=0.0d0
        return
    endif
    if(fcc.lt.faa) then
        r=dsqrt(fcc)
    else
        r=dsqrt(faa)
    endif
    dz=r/dfloat(nz)
    ichk=0
    do iz=1,nz+1
        z=dz*dfloat(iz-1)
        x=dsqrt(dabs(fcc-z**2))
        y=dsqrt((faa-z**2))
        term=x*z+y*z
        if(term.gt.fca)then
            ichk=1
            z2=z
            z1=z-dz
            exit
        endif
    enddo
    if(ichk.eq.0) then
        write(*,*)'no solution for c,alfa;det=',c,alfa,det
        return
    endif
    x1=dsqrt(fcc-z1**2)
    y1=dsqrt(faa-z1**2)
    x2=dsqrt(fcc-z2**2)
    y2=dsqrt(faa-z2**2)
    term1=x1*z1+y1*z1
    term2=x2*z2+y2*z2
    slope=(term2-term1)/(z2-z1)
    const=term1-slope*z1
    z=(fca-const)/slope
    x=dsqrt(fcc-z**2)
    y=dsqrt(faa-z**2)
    rcc=x
    rca=z
    raa=y

    return
end subroutine
!=============================================================================================
!    Another interpolation routine. When outside the grid, returns one of the boundary
!    points (not really the nearest endpoint, and only checks if greater than max inds).
!    Otherwise, again does bilinear interpolation.
!=============================================================================================
subroutine inttwo(c,alfa,quan,valu)
    use varbl
    implicit real*8 (a-h,o-z)
    dimension quan(icmax,iamax)
    ic=(c-ca(1,1))/dcc+1.00000001d0
    ia=(alfa-ca(1,2))/daa+1.00000001d0
    if(ic.gt.icmax-1.or.ia.gt.iamax-1)then
        if(ic.gt.icmax-1) then
            valu=quan(icmax,ia)
        endif
        if(ia.gt.iamax) then
            valu=quan(ic,iamax)
        endif
!csy  write(*,*)'error in inttwo','Q20=',c,'Q30=',alfa
        return
    endif

    p=(c-dcc*dfloat(ic-1)-qci)/dcc
    q=(alfa-daa*dfloat(ia-1)-qai)/daa
    f00=quan(ic,ia)
    f10=quan(ic+1,ia)
    f01=quan(ic,ia+1)
    f11=quan(ic+1,ia+1)
    valu=f00+p*(f10-f00)+q*(f01-f00)+p*q*(f00+f11-f01-f10)
    return
end subroutine

!=============================================================================================
!    When the outer turning line is supplied in 'prob.in', literally just populates
!    the array containing the starting points and their probabilities (technically,
!    must set the correct option in 'lan.in'). Otherwise, attempts to find the outer
!    turning line itself, which... just... why?
!=============================================================================================
subroutine set_langevin_starting_points(kotl,alfa1,alfa2,alfstp,cstart,ki)
    use varbl
    implicit real*8(a-h,o-z)

    if(kotl.eq.0)then
        allocate (cin(npro))
        allocate (ain(npro))
        allocate (pin(npro))
        do i=1,npro
            pin(i)=pro(i)
            cin(i)=cp(i)
            ain(i)=ap(i)
        enddo
        ki=npro
        return
    else
        allocate (cin(10000))
        allocate (ain(10000))
        allocate (pin(10000))
        nalf=(alfa2-alfa1)/alfstp+1
        alfa0=alfa1-2.0d0
        dalf0=0.1d0
        delc =0.2d0
        cdif =10.0d0
        cdff =1.0d0
        nc=(cc(icmax)-cstart)/0.2d0+1.0d0
        ki=0
        do ia=1,nalf
            alfa=alfa1+dfloat(ia-1)*alfstp
            do ic=1,nc
                c=cstart+dfloat(ic-1)*delc
                call inttwo(c,alfa,ff1,poten)
                if(poten.le.e0)then
                    if(ia.eq.1)cold=c
                    dif=dabs(cold-c)
                    if(dif.gt.cdif)then
                        ndf=dif/cdff
                        na=(aa(iamax)-alfa0)/dalf0+1.0d0
                        do ics=1,ndf
                            cii=c
                            if(c.lt.cold)cii=cold
                            cs=cii-dfloat(ics)*cdff
                            do ias=1,na
                                as=alfa0+dfloat(ias-1)*dalf0
                                as1=alfa0+dfloat(ias)*dalf0
                                call inttwo(cs,as,ff1,pots)
                                call inttwo(cs,as1,ff1,pots1)
                                if(pots.gt.e0.and.pots1.le.e0)then
                                    ki=ki+1
                                    cin(ki)=cs
                                    ain(ki)=as
                                    exit
                                endif
                            enddo
                        enddo
                    endif
                    cold=c
                    ki=ki+1
                    cin(ki)=c
                    ain(ki)=alfa
                    exit
                endif
            enddo
        enddo
        do i=1,ki
            ds=10.0
            pin(i)=1.0d0
            if(cin(i).gt.370.0)then
                pin(i)=0.0d0
                goto 219
            endif
            if(ain(i).gt.51.0d0*2.0.and.ain(i).lt.55.0d0*2.0)then
                pin(i)=pin(i-1)-0.010d0*(ain(i)-ain(i-1))
                goto 219
            elseif(ain(i).gt.55.0d0*2.0)then
                pin(i)=pin(i-1)*0.90d0
                goto 219
            endif
            do j=1,npro
                dd=sqrt((cin(i)-cp(j))**2+(ain(i)*0.5d0-ap(j))**2)
                if(dd.le.ds)then
                    ds=dd
                    pin(i)=pro(j)
                endif
            enddo
219         continue
            write(*,*)i,'of',ki,cin(i),ain(i)*0.5d0,pin(i)
        enddo
        return
    endif
end subroutine
!=============================================================================================
!       All subroutines below used *solely* to sample randomly from a normal distribution
!=============================================================================================
subroutine gran(idum,rand,x,y,dx,n)
!--- call subroutine table1 once before calling this routine
!--- generates one gaussian distributed random number=ran
! Uses inverse transform sampling, https://en.wikipedia.org/wiki/Inverse_transform_sampling
    implicit real*8 (a-h,o-z)
    external ran
    real ran
    dimension x(10001),y(10001)
    xran=ran(idum)
    call intpol2(xran,valu,x,y,dx,n)
    rand=valu
    return
end subroutine
!----------------------------------------------------------
subroutine intpol2(arg,valu,x,y,dx,n)
    implicit real*8 (a-h,o-z)
!-- for equispaced x-values at interval of dx at (n+1) points
!-- for a given value of x=arg, interpolates to find corresponding y=valu
!-- not for general purpose, special provision inserted at two end points
    dimension x(10001),y(10001)
    if(arg.lt.x(1).and.arg.ge.0.0d0) then
        valu=y(1)
        return
    endif
    if(arg.gt.x(n+1).and.arg.le.1.0d0)then
        valu=y(n+1)
        return
    endif
    i=(arg-x(1))/dx+1.00000001d0
    if(arg.le.x(i))then
        valu=y(i)
        return
    else
        if (arg.gt.x(i).and.arg.lt.x(i+1)) then
            slope=(y(i+1)-y(i))/(x(i+1)-x(i))
            valu=slope*(arg-x(i))+y(i)
            return
        else
            write(*,*)arg,'argument outside interpolation range in intpol2'
        endif
    endif
    return
end subroutine
!=============================================================================================
!    A random number generator taken from e.g.
!    https://en.wikipedia.org/wiki/Combined_linear_congruential_generator#Example
!    or
!    https://www.mi.fu-berlin.de/inf/groups/ag-tech/teaching/2012_SS/L_19540_Modeling_and_Performance_Analysis_with_Simulation/06.pdf (slide 31)
!=============================================================================================
function ran(idum)
    integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real ran,AM,EPS,RNMX
    parameter IM1=2147483563,IM2=2147483399,AM=1./IM1,&
        IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS
    integer idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if(idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if (iy.lt.1)iy=iy+IMM1
    ran=min(AM*iy,RNMX)
    return
end function
!------------------------------------------------
subroutine table1(b,x,y,dx,n)
!--- calculates the (x,y) table at n+1 points for a gaussian
!---- x-values are equispaced at interval of dx
!--- sqrt(b/pi)*exp(-b*y*y)
    implicit real*8 (a-h,o-z)
    dimension g(2001)
    dimension x(n+1),y(n+1)
    dimension xx(n+1),yy(n+1)

    pi=4.0D0*datan(1.0D0)
    m=2000

!--- n=no. of points at which F(y)(=x) is calculated as a function of y
!--- m=no. of integration points to calculate F(y)
    const=dsqrt(b/pi)
    ymin=-dsqrt(15.0d0/b)
    ymax= dsqrt(15.0d0/b)
    dy=(ymax-ymin)/dfloat(n)
    do i=1,n+1
        yy(i)=ymin+dfloat(i-1)*dy
        y1=0.0d0
        y2=yy(i)
        h=(y2-y1)/dfloat(m)
        do j=1,m+1
            yrun=y1+dfloat(j-1)*h
            g(j)=const*dexp(-b*yrun*yrun)
        enddo
        call simp(g,h,m,valu)
        xx(i)=0.5d0+valu
    enddo

!-- the table (xx(i),yy(i)) is now ready
!-- another table (x,y) will be now prepared where x-values are equispaced
    x1=xx(1)
    x2=xx(n+1)
    dx=(x2-x1)/dfloat(n)
    do i=1,n+1
        arg=x1+dfloat(i-1)*dx
        call intpol1(arg,valu,xx,yy,n)
        x(i)=arg
        y(i)=valu
    enddo
    return
end subroutine

!-----------------------------------------------------------
subroutine intpol1(arg,valu,x,y,n)
    implicit real*8 (a-h,o-z)
!-- for non-equispaced x-values at (n+1) points
!-- for a given value of x=arg, interpolates to find corresponding y=valu
!-- not for general purpose, special provision inserted at two end points
    dimension x(n+1),y(n+1)
    if(arg.lt.x(1).and.arg.ge.0.0d0) then
        valu=y(1)
        return
    endif
    if(arg.gt.x(n+1).and.arg.le.1.0d0)then
        valu=y(n+1)
        return
    endif
    i=1
3   if(arg.eq.x(i)) then
        valu=y(i)
        return
    else if (arg.gt.x(i).and.arg.lt.x(i+1)) then
        slope=(y(i+1)-y(i))/(x(i+1)-x(i))
        valu=slope*(arg-x(i))+y(i)
        return
    else
        i=i+1
        m=n+1
        if(i.eq.m.and.arg.eq.x(m)) go to 6
        if(i.gt.n.and.i.ne.m) go to 5
        go to 3
    endif
6   valu=y(m)
    return
5   write(*,*)arg,'argument outside interpolation range in intpol1'
    return
end subroutine

!----------------------------------------------------------
subroutine simp(g,h,n,valu)
    implicit real*8 (a-h,o-z)
!INTEGRATE FUNCTION g(z) using (n+1) equispaced points with step h
    DIMENSION g(n+1)

    SUM1=0.d0
    sum2=0.d0
    if(n.gt.2)then
        do i=2,n,2
            sum1=sum1+g(i)
        enddo
        do j=3,n-1,2
            sum2=sum2+g(j)
        enddo
    endif
    valu=(h/3.d0)*(g(1)+4.d0*sum1+2.d0*sum2+g(n+1))
    RETURN
END subroutine
