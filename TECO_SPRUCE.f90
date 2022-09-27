    program TECO_MCMC
        
! USE IFPORT !! Was: USE DFPORT
!   for parameter file
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,Rl0,Rs0,Rr0
    real paraest(19,40000)
    integer seq,Pselect
    character(len=120) parafile,daparfile,outdir
    character(len=150) paraestfile
    integer,parameter :: partotal=35
    integer,dimension(partotal):: DApar
    real,dimension(partotal) :: parval,parmin,parmax
    character(len=250) indexstring
    
!   for climate file
    integer, parameter :: ilines=150000
    integer, parameter :: iiterms=7
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real forcing_data(iiterms,ilines)
    character(len=150) climatefile
    integer yr_length,lines
    integer,dimension(ilines):: year_seq1,doy_seq1,hour_seq1
    real forcing_data1(iiterms,ilines)
    character(len=150) climatefile1
    integer yr_length1,lines1    
    integer,dimension(ilines):: year_seq2,doy_seq2,hour_seq2
    character(len=150) climatefile2,forcingdir
    real forcing_data2(iiterms,ilines)
    integer yr_length2,lines2
    
!   for observation files    
    real obs_spruce(12,1000),std(12,1000)
    real treatment
    character(len=150) obsfile1,obsfile2,covfile
    integer len1,len2
    
    
!   for MCMC
    integer MCMC ! 0:run model simulation; 1:data assimilation 
    integer IDUM,upgraded,isimu
!    integer, parameter :: npara=18       ! Number of parameters to be estimated
    integer npara
    real search_length
    real J_last
    real Simu_dailyflux(12,10000)
!    real coef(npara),coefac(npara),coefnorm(npara).
!    real coefmax(npara),coefmin(npara)
    real r,fact_rejet
!    real gamma(npara,npara),gamnew(npara,npara)     ! covariance matrix
    real, allocatable :: coef(:), coefac(:), coefnorm(:)
    real, allocatable :: coefmax(:),coefmin(:)
    real, allocatable :: gamma(:,:),gamnew(:,:)
    integer,allocatable :: coefindex(:)
    
    integer k1,k2,rejet,paraflag,k3
    integer, parameter :: nc=100
    integer, parameter :: ncov=500
!! nc: the multiplicative constant will be adjusted !!
!!        every nc iterations, to preserve an adequate!!
!!        acceptation rate   
!! ncov: the covariance matrix gamma will be updated  !!
!!        every ncov iterations		
!    real coefhistory(ncov,npara)
    real, allocatable :: coefhistory(:,:)
    character(len=150) outfile,MCMCargu,yrargu,dyargu
    character(len=150) Targu,CO2argu
    
!   for consts parameteres
    real,dimension(3):: tauL,rhoL,rhoS
    real pi,emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk
    real H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta
    real conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm
    real Entrpy,gam0,gam1,gam2
    
!   for initialize
    real fwsoil,topfws,omega,Storage,nsc
    real fwsoil_initial,topfws_initial,omega_initial
    real Storage_initial,nsc_initial
    real wcl(10),QC(8)
    real wcl_initial(10),QC_initial(8)
    integer yrs_eq,rep,yrlim,dylim
    character(len=150) my_fmt
    
 
    yrlim = 2014
    dylim = 365
    Ttreat = 0.0
    CO2treat = 380.0
!   Read parameters from file
    call getarg(1,parafile)
    !parafile='input/SPRUCE_pars.txt'
    call Getparameters(lat,longi,wsmax,wsmin,           &              
    &   LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,           &
    &   SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n,         &
    &   a1,Ds0,Vcmax0,extkU,xfang,alpha,                 &
    &   Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,         &
    &   Tau_Micro,Tau_slowSOM,Tau_Passive,              &
    &   gddonset,Q10,RL0,Rs0,Rr0,parafile)
    parval = (/lat,longi,wsmax,wsmin,           &              
    &   LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,           &
    &   SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n,         &
    &   a1,Ds0,Vcmax0,extkU,xfang,alpha,                 &
    &   Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,         &
    &   Tau_Micro,Tau_slowSOM,Tau_Passive,              &
    &   gddonset,Q10,RL0,Rs0,Rr0/)
    
    
!   Read climatic forcing
!    climatefile='SPRUCE_forcing.txt'
    call getarg(2,climatefile1)
    !climatefile1='input/SPRUCE_forcing.txt'
    call Getclimate(year_seq1,doy_seq1,hour_seq1,          &
    &   forcing_data1,climatefile1,lines1,yr_length1)

!   Read observation data
    call getarg(3,obsfile1)
    !obsfile1='input/SPRUCE_obs.txt'
    treatment=0.    ! Ambient temperature
    call GetObsData(obs_spruce,std,len1,obsfile1)      
            
!   initiations for canopy model, including canopy traits variation in a year
    call consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,Vcmax0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)
    fwsoil=1.0
    topfws=1.0
    omega=1.0
    do i=1,10
        wcl(i)=wsmax/100.
    enddo 
    Storage=32.09           !g C/m2
    nsc=85.35
!         put the values into a matrix
!   QC=(/300.,6300.,300.,119.,300.,322.,8834.,312./)
!   QC=(/100.,800.,100.,39.,100.,122.,834.,12./)
 !   QC=(/440.,700.,300.,119.,300.,322.,38340.,23120./)
    QC=(/300.,650.,100.,119.,300.,322.,38340.,23120./)
    
    
!   Start main loop
    call getarg(4,outdir)
    !outdir = 'output'
    write(outfile,"(A120,A18)") trim(outdir),"/SPRUCE_yearly.txt"
    outfile = trim(outfile)
    outfile = adjustl(outfile)
    open(61,file=outfile)
    write(outfile,"(A120,A22)") trim(outdir),"/Simu_dailyflux001.txt"
    outfile = trim(outfile)
    outfile = adjustl(outfile)
    open(62,file=outfile)
    
    call getarg(5,MCMCargu)
    read(MCMCargu,'(i1)') MCMC
    !MCMC = 2

    call getarg(6,DAparfile)
    !DAparfile='input/SPRUCE_da_pars.txt'
    call GetDAcheckbox(DApar,parmin,parmax,DAparfile)

    if(MCMC.eq.1) GOTO 100
    if(MCMC.eq.2) GOTO 150

    year_seq = year_seq1
    doy_seq = doy_seq1
    hour_seq = hour_seq1
    forcing_data = forcing_data1
    climatefile = climatefile1
    lines = lines1
    yr_length = yr_length1
    yrs_eq=yr_length*0  ! spin up length 
    call TECO_simu(MCMC,Simu_dailyflux,      &
     &        obs_spruce,yrlim,dylim,Ttreat,CO2treat,              &
     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,    &
     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
     &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
     &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
     &        tau_SlowSOM,tau_Passive,gddonset,                     &
     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)    
    write(*,*)'run simulation'
    return
150 continue    
    !   Update posterior parameters
    write(paraestfile,"(A120,A12)") trim(outdir),"/Paraest.txt"
    paraestfile = trim(paraestfile)
    paraestfile = adjustl(paraestfile)
    call Getparaest(paraestfile,paraest,seq,npara,indexstring)
    allocate(coefindex(npara))
    write (my_fmt, '(a,i0,a)') '(',npara,'I12)'
    read(indexstring,my_fmt) coefindex

    call getarg(8,yrargu)
    read(yrargu,'(i4)') yrlim
    !yrlim = 2024
    call getarg(9,dyargu)
    read(dyargu,'(i3)') dylim
    !dylim = 365
    call getarg(10,Targu)
    read(Targu,'(f9.3)') Ttreat
    !Ttreat = 0.0
    call getarg(11,CO2argu) 
    read(CO2argu,'(f9.3)') CO2treat
    !CO2treat = 380.0
    
    
    DO rep=1,100
        
    CALL random_number(randnum)
    Pselect = int(seq/2+randnum*(seq-seq/2))
    
    !Pselect = 10000
    do k1=1,npara
        parval(coefindex(k1))=paraest(k1+1,Pselect)
    enddo
        SLA = parval(12)
        GLmax = parval(13)
        GRmax = parval(14)
        Gsmax = parval(15)
        Vcmax0 = parval(19)
        Tau_Leaf = parval(23)
        Tau_Wood = parval(24)
        Tau_Root = parval(25)
        Tau_F = parval(26)
        Tau_C = parval(27)
        Tau_Micro = parval(28)
        Tau_slowSOM = parval(29)
        Tau_Passive = parval(30)
        gddonset = parval(31)
        Q10 = parval(32)
        RL0 = parval(33)
        Rs0 = parval(34)
        Rr0 = parval(35)
    
!   Read generated climatic forcing
    call getarg(7,forcingdir)
!    forcingdir = 'input/Weathergenerate'
    write(climatefile2,"(A120,A10,I3.3,A4)") trim(forcingdir),"/EMforcing",rep,".csv"
    climatefile2=trim(climatefile2)
    climatefile2=adjustl(climatefile2)
    call Getclimate(year_seq,doy_seq,hour_seq,          &
    &   forcing_data,climatefile2,lines,yr_length)
    do k1=1,lines1
        year_seq(k1)=year_seq1(k1)
        doy_seq(k1)=doy_seq1(k1)
        hour_seq(k1)=hour_seq1(k1)
        forcing_data(1,k1)=forcing_data1(1,k1)
        forcing_data(2,k1)=forcing_data1(2,k1)
        forcing_data(3,k1)=forcing_data1(3,k1)
        forcing_data(4,k1)=forcing_data1(4,k1)
        forcing_data(5,k1)=forcing_data1(5,k1)
        forcing_data(6,k1)=forcing_data1(6,k1)
        forcing_data(7,k1)=forcing_data1(7,k1)
    enddo
    
    write(outfile,"(A120,A15,I3.3,A4)") trim(outdir), "/Simu_dailyflux",rep,".txt"
    outfile=trim(outfile)
    outfile=adjustl(outfile)
    open(62,file=outfile)
    fwsoil=1.0
    topfws=1.0
    omega=1.0
    do i=1,10
        wcl(i)=wsmax/100.
    enddo 
    Storage=32.09           !g C/m2
    nsc=85.35
!         put the values into a matrix
!   QC=(/300.,6300.,300.,119.,300.,322.,8834.,312./)
!   QC=(/100.,800.,100.,39.,100.,122.,834.,12./)
 !   QC=(/440.,700.,300.,119.,300.,322.,38340.,23120./)
    QC=(/300.,650.,100.,119.,300.,322.,38340.,23120./)
    
    yrs_eq=yr_length*0  ! spin up length 
    call TECO_simu(MCMC,Simu_dailyflux,      &
     &        obs_spruce,yrlim,dylim,Ttreat,CO2treat,              &
     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,    &
     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
     &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
     &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
     &        tau_SlowSOM,tau_Passive,gddonset,                     &
     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)    
    write(*,*)'run forecasting',rep
    close(62)
    enddo ! END of rep

    return
    
100 continue

    write(outfile,"(A120,A12)") trim(outdir),"/Paraest.txt"
    outfile = trim(outfile)
    outfile = adjustl(outfile)
    open(71,file=outfile)
    year_seq = year_seq1
    doy_seq = doy_seq1
    hour_seq = hour_seq1
    forcing_data = forcing_data1
    climatefile = climatefile1
    lines = lines1
    yr_length = yr_length1
    
    fwsoil_initial = fwsoil
    topfws_initial = topfws
    omega_initial = omega
    wcl_initial = wcl
    Storage_initial = Storage
    nsc_initial = nsc
    QC_initial = QC
!   end of spin up
    
    npara=sum(DApar)
    allocate(coef(npara),coefac(npara),coefnorm(npara))
    allocate(coefindex(npara))
    allocate(coefmax(npara),coefmin(npara))
    allocate(gamma(npara,npara),gamnew(npara,npara))
    allocate(coefhistory(ncov,npara))
!	simulation begins here
    J_last=9000000.0
    IDUM = 542
    upgraded=0
    new=0
    k3=0
    j=0
    do i=1,35
       if (DApar(i).eq. 1) then
           j=j+1
           coef(j)=parval(i)
           coefindex(j)=i
           coefmin(j)=parmin(i)
           coefmax(j)=parmax(i)
       endif
    enddo
    write(71,*) npara
    write(71,*)(coefindex(i),i=1,npara)
    ! initialize covariance matrix
    covexist=0
    if(covexist.eq.1)then      ! If prior covariance exists, read from file
   
        write(covfile,"(A120,A15)") trim(outdir),"/covariance.txt"
        covfile = trim(covfile)
        covfile = adjustl(covfile)
        call getCov(gamma,covfile,npara)

        call racine_mat(gamma,gamnew,npara)      ! square root of covariance matrix
        gamma=gamnew
        do k1=1,npara
!           coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
            coefnorm(k1)=0.5
            coefac(k1)=coefnorm(k1)
        enddo
    else
        coefac=coef
    endif
 
    fact_rejet=2.4/sqrt(real(npara))
    search_length=0.05
    rejet = 0
    

    do isimu=1,50000
!       generate parameters
        if(covexist.eq.1)then 
            paraflag=1
            do while(paraflag.gt.0)
                call gengaussvect(fact_rejet*gamma,coefac,coefnorm,npara)
                paraflag=0
                do k1=1,npara
                    if(coefnorm(k1).lt.0. .or. coefnorm(k1).gt.1.)then
                    paraflag=paraflag+1
                    write(*,*)'out of range',paraflag
                    endif
                enddo
            enddo
            do k1=1,npara
                coef(k1)=coefmin(k1)+coefnorm(k1)*(coefmax(k1)-coefmin(k1))
            enddo
        else
            call coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        endif

!         update parameters
        do k1=1,npara
            parval(coefindex(k1))=coef(k1)
        enddo
        SLA = parval(12)
        GLmax = parval(13)
        GRmax = parval(14)
        Gsmax = parval(15)
        Vcmax0 = parval(19)
        Tau_Leaf = parval(23)
        Tau_Wood = parval(24)
        Tau_Root = parval(25)
        Tau_F = parval(26)
        Tau_C = parval(27)
        Tau_Micro = parval(28)
        Tau_slowSOM = parval(29)
        Tau_Passive = parval(30)
        gddonset = parval(31)
        Q10 = parval(32)
        RL0 = parval(33)
        Rs0 = parval(34)
        Rr0 = parval(35)

        yrs_eq = 0
        call TECO_simu(MCMC,Simu_dailyflux,            &
        &        obs_spruce,yrlim,dylim,Ttreat,CO2treat,              &
        &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
        &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,         &
        &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
        &        Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,      &
        &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
        &        tau_Leaf,tau_Wood,tau_Root,tau_F,tau_C,tau_Micro,     &   ! the unit is year
        &        tau_SlowSOM,tau_Passive,gddonset,                     &
        &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)   

!         M-H algorithm
        tmp_up=upgraded
        
        call costFObsNee(Simu_dailyflux,   &
            &   obs_spruce,std,len1,   &
            &   J_last,upgraded)

        if(upgraded.gt.tmp_up)then 
            new=new+1
            if(covexist.eq.1)then
                coefac=coefnorm
                coefhistory(new,:)=coefnorm
            else
                coefac=coef
                do k1=1,npara
                    coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
                enddo
            endif
            coefhistory(new,:)=coefnorm
            if(new.ge.ncov)new=0
            write (my_fmt, '(a,i0,a)') '(i12,",",',npara,'(F15.4,","))'
            write(71,my_fmt)upgraded,(coef(i),i=1,npara)
            if(upgraded.gt.1000 .and. k3.lt.500)then
                CALL random_number(r)
                if(r.gt.0.95)then
                    k3=k3+1
                    write(outfile,"(A120,A15,I3.3,A4)") trim(outdir), "/Simu_dailyflux",k3,".txt"
                    outfile=trim(outfile)
                    outfile=adjustl(outfile)
                    open(62,file=outfile)
                    do i=1,1460
                        write(62,602)i,(Simu_dailyflux(j,i),j=1,12)
                    enddo
602                 format((i7),",",12(f15.4,","))
                    close(62)
                    
                endif
            endif
        else
            reject=reject+1
        endif

        fwsoil = fwsoil_initial
        topfws = topfws_initial
        omega = omega_initial
        wcl = wcl_initial
        Storage=Storage_initial
        nsc=nsc_initial
        QC = QC_initial
      
        write(*,*)isimu,upgraded
    
    	! updates of the multiplicative constant
        if(covexist.eq.1)then
            if(mod(isimu,nc).eq.0)then
                if ((1. - real(rejet)/real(nc)) < 0.23) then
                !    fact_rejet = fact_rejet*0.9
                else
                    if ((1. - real(rejet)/real(nc)) > 0.44) then
                !    fact_rejet = fact_rejet * 1.1
                    endif
                endif
            rejet=0
            write(*,*)'search length is', search_length
            endif
        else
            if(mod(isimu,nc).eq.0)then
                if(real(upgraded)/real(isimu) .lt. 0.23)then
                !    search_length=search_length*0.9
                else
                    if(real(upgraded)/real(isimu) .gt. 0.44)then
                !        search_length=search_length*1.1
                    endif
                endif
                reject=0
                write(*,*)'search length is', search_length
            endif
        endif

	! updates of the covariance matrix
        if(covexist.eq.0 .and. mod(upgraded,ncov).eq.0)then
            covexist=1
            coefac=coefnorm
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew
                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
        endif
	if (mod(upgraded,ncov).eq.0 .and. covexist.eq.1) then
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew
                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
	endif

    enddo !isimu
    !write(outfile,"(A120,A21)") trim(outdir),"/covvariance_temp.txt"
    !outfile = trim(outfile)
    !outfile = adjustl(outfile)
    !open(72,file=outfile)
    !do i=1,npara
    !    write(72,*) (gamma(j,i),j=1,npara)
    !enddo
    !close(72)    
    close(61)

    close(71)
    
    deallocate(coef,coefac,coefnorm)
    deallocate(coefindex)
    deallocate(coefmax,coefmin)
    deallocate(gamma,gamnew)
    deallocate(coefhistory)
    
    write(*,*)'run MCMC'
    end


! ====================================================================
    subroutine TECO_simu(MCMC,Simu_dailyflux,      &
     &        obs_spruce,yrlim,dylim,Ttreat,CO2treat,        &
     &        forcing_data,yr_length,year_seq,doy_seq,hour_seq,lines,   &
     &        fwsoil,topfws,omega,wcl,Storage,nsc,yrs_eq,QC,            &
     &        lat,longi,wsmax,wsmin,LAIMAX,LAIMIN,rdepth,     &
     &        Rootmax,Stemmax,SapR,SapS,SLAx,GLmx,GRmx,Gsmx,      &
     &        stom_n,a1,Ds0,Vcmax0,extkU,xfang,alpha,               &
     &        tau_L,tau_W,tau_R,tau_F,tau_C,tau_Micr,     &   ! the unit is year
     &        tau_Slow,tau_Pass,gddonset,                     &
     &        Q10,Rl0,Rs0,Rr0,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &    Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,&
     &    wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,&
     &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)
      
      implicit none
      
   !==================== test variables
    real fwsoil_yr,omega_yr,topfws_yr, difference,diff_yr,diff_d
            
      integer, parameter :: iiterms=7            ! 9 for Duke forest FACE
      integer, parameter :: ilines=150000         ! the maxmum records of Duke Face, 1998~2007
      real, parameter:: times_storage_use=720.   ! 720 hours, 30 days
      integer  lines,idays,MCMC
      integer,dimension(ilines):: year_seq,doy_seq,hour_seq
      real forcing_data(iiterms,ilines),input_data(iiterms,ilines)
      real Simu_dailyflux(12,10000)
      real obs_spruce(12,1000)
      integer pheno,phenoset
!      site specific parameters
      real lat,longi,rdepth,LAIMAX,LAIMIN
      real wsmax,wsmin,co2ca,CO2treat
      real tau_L,tau_W,tau_R
      real tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass
      real TauC(8)
!      the variables that should be initialized in the begining
      real Q_soil
      real QC(8) !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
      real Pool1,Pool2,Pool3,Pool4,Pool5,Pool6,Pool7,Pool8
      real out1_yr,out2_yr,out3_yr,out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
      real OutC(8)
      real Rh_pools(5)
!      for soil conditions
      real WILTPT,FILDCP,infilt
      real Rsoilabs
      real fwsoil,topfws,omega
!      for plant growth and allocation
      real NSC,NSCmin,NSCmax,add               ! none structural carbon pool
      real Growth,Groot,Gshoot,GRmax           ! growth rate of plant,root,shoot,and max of root
      real St,Sw,Ss,Sn,Srs,Sps,fnsc,Weight     ! scaling factors for growth
!      variables for canopy model
      real evap,transp,ET,G
      real wind,eairp,esat,rnet
      real Pa_air 
      real gpp,gpp_ra,NPP,NEE,NEP,gpp_d,NPP_d
      real evap_d,transp_d
      real,dimension(3):: tauL,rhoL,rhoS,reffbm,reffdf,extkbm,extkdm
      real,dimension(2):: Radabv
      real Qcan(3,2)
!      parameters for photosynthesis model
      real stom_n,a1,Ds0,Vcmx0,Vcmax0,extkU,xfang,alpha
      real pi,emleaf,emsoil
      real Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat
      real wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci
      real Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2
!     for nitrogen sub-model
      real CNmin,CNmax,NSNmax,NSNmin
      real NSN
!      QNleaf,QNwood,QNroot,QNfine,QNcoarse,QNmicr,QNslow,QNpass
!      CN_leaf,CN_wood,CN_root,CN_fine,CN_coarse,CN_micr,CN_slowC,CN_pass
      real QN(8),CN0(8),CN(8),OutN(8),QNplant,QNminer
      real N_leaf,N_wood,N_root,N_deficit
      real N_LF,N_WF,N_RF
      real N_uptake,N_leach,N_vol,N_fixation,N_deposit,N_fert
      real N_up_d,N_fix_d,N_dep_d,N_leach_d,N_vol_d
      real N_up_yr,N_fix_yr,N_dep_yr,N_leach_yr,N_vol_yr
      real N_miner,alphaN
      real SNvcmax,SNgrowth,SNRauto,SNrs

!      additional arrays to allow output of info for each layer
      real,dimension(5):: RnStL,QcanL,RcanL,AcanL,EcanL,HcanL
      real,dimension(5):: GbwcL,GswcL,hG,hIL
      real,dimension(5):: Gaussx,Gaussw,Gaussw_cum 
!      for phenology
      real LAI,bmroot,bmstem,bmleaf,bmplant,totlivbiom,ht
      real SLA,SLAx,L_fall,L_add,litter,seeds
      real GDDonset,GDD5,accumulation,storage,stor_use,store
      real RaL,RaS,RaR  !allocation to respiration
      real alpha_L,alpha_W,alpha_R ! allocation ratio to Leaf, stem, and Root
      real Q10,Rl0,Rs0,Rr0         ! parameters for auto respiration
      real Rgrowth,Rnitrogen,Rmain,Rauto !respirations
      real RmLeaf,RmStem,RmRoot          ! maintanence respiration
      real RgLeaf,RgStem,RgRoot          ! growth respiration
      real RaLeaf,RaStem,RaRoot
      real Rsoil,Rhetero,Rtotal
      real Ra_Nfix,Rh_Nfix
      real gpp_yr,NPP_yr,NEE_yr,RaL_yr,RaR_yr,RaS_yr,Rh_yr
      real Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,Ra_yr
      real R_Ntr_yr
      real NPPL_yr,NPPR_yr,NPPS_yr,NPP_L,NPP_R,NPP_W
      real Rootmax,Stemmax,SapS,SapR,StemSap,RootSap
      REAL ws,wdepth
!      climate variables for every day
      real Ta,Tair,Ts,Tsoil,Ttreat
      real doy,hour,Dair,Rh,radsol
      real PAR
!      output daily means of driving variables
      real CO2air_d_avg,SWdown_d_avg,Psurf_d_avg
      real Rain_d_avg,Tair_d_avg,Wind_d_avg

!      output from canopy model
      real evap_yr,transp_yr
      real,dimension(10):: thksl,wupl,evapl,wcl,FRLEN
      real runoff,runoff_d,runoff_yr,rain,rain_d,rain_yr
      real wsc(10),ws1,ws2,dws,net_dws
      real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop
      real Hcanop,Hcanop_d
      real Raplant,Glmax,Gsmax,Rh_d
      real GLmx,Gsmx,GRmx
!      output for ORNL model comparison
      real CO2h,PARh,ATh,STh,VPDh,SWh
      real RECOh
      real ETh,Th,Eh,INTh,ROh,DRAINh,LEh,SHh
      real LWH,Rgrowth_d,abvLitter,blLitter
!     daily output
      real PAR_d,AT_d,ST_d,VPD_d
      real SW_d,NEP_d,NEE_d,RECO_d
      real Ra_d,RLEAV_d,RWOOD_d,RROOT_d,RHET_d,RSOIL_d,ET_d,T_d
      real E_d,INT_d,RO_d,DRAIN_d,LE_d,SH_d,CL_d,CW_d,CFR_d,TNC_d
      real CSOIL_d,GL_d,GW_d,GR_d,LFALL_d,LMA_d,NCAN_d,NWOOD_d
      real GL_yr,GR_yr,GW_yr
      real NFR_d,NSOIL_d,NUP_d,NMIN_d,NVOL_d,NLEACH_d
      real N_LG_d,N_WG_d,N_RG_d
      real N_LF_d,N_WF_d,N_RF_d
      real WFALL_D,RFALL_D
      real Simu_lit

!      NEE observation
      real NEE_annual,Cumol2gram
      real NEE_annual_array(30)
      integer year_array(30),year_obs
!     for loops
      integer jrain,W_flag(7)
      integer onset !flag of phenological stage
      integer year,yr,days,i,j,m,n,yrs_eq,hoy,iyr,daily
      integer k1
      integer lines_NEE,yr_NEE
      integer istat1,istat2,istat3,istat4
      integer dtimes,yr_length
      integer num_scen,isite
      integer idoy,ihour,ileaf,first_year
      integer dylim,yrlim


!     Default C/N ratios of Oak Ridge FACE
!      data CN0 /45.,350.,60.,40.,300.,10.,15.,8./
!     Default C/N ratios of Duke FACE
      data CN0 /50.,350.,60.,40.,300.,10.,20.,12./
!     thickness of every soil layer
      data thksl /10.,10.,10.,10.,10.,20.,20.,20.,20.,20./
!     ratio of roots in every layer, Oak Ridge FACE
      data FRLEN /0.1,0.25,0.25,0.2,0.1,0.05,0.025,0.015,0.005,0.005/

!     Nitrogen input
!      N_deposit=0.000144634702 !(gN/h/m2, 1.2+0.067 gN/yr/m2,Oak ridge)
!     0.7 gN/yr/m2, 13.4 kg N ha-1 yr-1, 2000, Dentener et al. 2006, GBC, Duke FACE
      N_deposit=2.34/8760. !(gN/h/m2, )

!      N_fert=0. ! (20.0 gN m-2 yr-1, in spring, from 2004, Oak Ridge)
      N_fert=0. !5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)

!         the unit of residence time is transformed from yearly to hourly
          tauC=(/tau_L,tau_W,tau_R,tau_F,tau_C,&
     &           tau_Micr,tau_Slow,tau_Pass/)*8760.
     
          SLA=SLAx/10000.         ! Convert unit from cm2/g to m2/g
!         growth rates of plant
          GLmax=GLmx/8760.
          GRmax=GRmx/8760.
          Gsmax=GSmx/8760.
!         end of setting parameters

 
          input_data=forcing_data
!         end of reading forcing data

!     ===============================================================
!         cycle  ! skip the following blocks, read input data only.
!         ===================================================
!         Initialize parameters and initial state:


          WILTPT=wsmin/100.0
          FILDCP=wsmax/100.0
!         define soil for export variables for satisfying usage of canopy submodel first time
!          wscontent=WILTPT

          infilt=0.

!         gddonset=320.0
          
          stor_use=Storage/times_storage_use
          accumulation=0.0
          SNvcmax=1.0


          LAI=LAIMIN
          bmleaf=QC(1)/0.48
          bmstem=QC(2)/0.48
          bmroot=QC(3)/0.48
          bmplant=bmstem+bmroot+bmleaf

!         initial values of Nitrogen pools and C/N ratio
          alphaN=0.0    ! the transfer of N before littering

          NSN=6.0
          QNminer= 1.2
          N_deficit=0
          CN=CN0
          QN=QC/CN0
          QNplant  =QN(1) + QN(2) + QN(3)

!=============================================================
          m=1
          n=1
          k1=1
          iyr=0
          idays=365
          daily=0
          first_year=2011
          !yr_length=1
          do yr=1,yrs_eq+yr_length  ! how many years
              if(yr.gt.3)then
!                  write(*,*)'One year done'
              endif
!            using ambient data to run equilibiurm, elevated only for the last cycle
             iyr=iyr+1
             if(iyr>yr_length)iyr=1

!!          leap year
           if(MOD(first_year+iyr-1,4).eq.0)then
                 idays=366
           else
                 idays=365
           endif
           
           
             GDD5=0.0
             onset=0
             phenoset=0
             diff_yr=0.0
             gpp_yr=0.0
             R_Ntr_yr=0.
             NPP_yr=0.0
             Rh_yr =0.0
             Rh4_yr=0.0
             Rh5_yr=0.0
             Rh6_yr=0.0
             Rh7_yr=0.0
             Rh8_yr=0.0
             Ra_yr =0.0
             GL_yr=0.0
             GW_yr=0.0
             GR_yr=0.0
             Pool1=0.0
             Pool2=0.0
             Pool3=0.0
             Pool4=0.0
             Pool5=0.0
             Pool6=0.0
             Pool7=0.0
             Pool8=0.0
             out1_yr=0.0
             out2_yr=0.0
             out3_yr=0.0
             out4_yr=0.0
             out5_yr=0.0
             out6_yr=0.0
             out7_yr=0.0
             out8_yr=0.0             
             NEE_yr=0.0
!            water fluxes
             rain_yr=0.0
             transp_yr=0.0
             evap_yr=0.0
             runoff_yr=0.0
             Simu_lit=0.
             
!            Nitrogen fluxes
             N_up_yr=0
             N_fix_yr=0.
             N_dep_yr=0.
             N_leach_yr=0.
             N_vol_yr=0.
             
        !============================== test variable
             fwsoil_yr=0.
             omega_yr=0.
             topfws_yr=0.
             
             hoy=0

!!     end of leap year
             do days=1,idays !the days of a year

!             Nitrogen fertilization since 2004 in Oak Ridge
!              if(yr>yrs_eq+5.and.days==135)then
!                  QNminer=QNminer+N_fert     !(20 gN/yr/m2,N fertiliztion in Spring)
!              endif

!             Nitrogen fertilization since 1999 in Duke
              if(yr>yrs_eq+1.and.(days==75.OR.days==105))then
                  QNminer=QNminer+N_fert     !(5.6 gN/yr/m2,N fertiliztion in March and Apr)
              endif

              StemSap=AMIN1(Stemmax,SapS*bmStem)   ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
              RootSap=AMIN1(Rootmax,SapR*bmRoot)
              NSCmin=5. 
              NSCmax=0.05*(StemSap+RootSap+QC(1))
              if(Ta.gt.5.0)GDD5=GDD5+Ta
!             THE FIRST PART:  coupled canopy and soil model
              diff_d = 0.0
              gpp_d   =0.0   ! daily
              gpp_ra  =0.0   ! daily
              NPP_d   =0.0   ! daily
              NEP_d=0.0
              NEE_d = 0.0
!             rain_d,transp_d,evap_d
              transp_d=0.0   ! daily
              Hcanop_d=0.0   ! daily
              evap_d  =0.0   ! daily
              ta=0.0         ! daily 
              Ts=0.0         ! daily
              rain_d=0.0     ! daily
              runoff_d=0.0    ! daily
              LE_d=0.0
              RaL=0.0
              RaS=0.0
              RaR=0.0
              Rauto=0.0
              Rh_d=0.0
              N_up_d=0.
              N_fix_d=0.
              N_dep_d=0.
              N_leach_d=0.
              N_vol_d=0.
              PAR_d=0.
              VPD_d=0.0
              RECO_d=0.0
              RLEAV_d=0.0 
              RWOOD_d=0.0
              RROOT_d=0.0
              GL_d   =0.0
              GW_d   =0.0
              GR_d   =0.0
              LFALL_d=0.0
              NUP_d=0.0
              NVOL_d=0.
              NLEACH_d=0.0
              NMIN_d=0.0
              N_LG_d=0.0
              N_WG_d=0.0
              N_RG_d=0.0
              N_LF_d=0.0
              N_WF_d=0.0
              N_RF_d=0.0
              WFALL_d=0.0
              RFALL_d=0.0
              dtimes=24 !how many times a day,24 means every hour
              do i=1,dtimes
                  if(i.gt.11)then
!                      write(*,*)'pause'
                  endif
                    
!                 input data
                  if(m > lines)then    ! Repeat forcing data for the whole time period
                !  if(yr.ge.1 .and. mod(days,idays).eq.0     &
                !  &     .and. mod(i,dtimes).eq.0)then   ! Repeat forcing data for specific time period
                     m=1		! m is row sequence in the file of input_data
                     n=1
                     hoy=0
                  endif
                  
                  year =year_seq(m)
                  doy  =doy_seq(m)
                  hour =hour_seq(m)+1
                  
                  !!       for Duke Forest
                  Tair=input_data(1,m)   ! Tair
                  Tsoil=input_data(2,m)    ! SLT
                  co2ca=380.0*1.0E-6
                  if (yr .gt. 4)then
                      Tair = Tair + Ttreat
                      Tsoil = Tsoil + Ttreat
                      co2ca=CO2treat*1.0E-6 ! CO2 concentration,ppm-->1.0E-6
                  endif                  
                  RH=input_data(3,m)
                  Dair=input_data(4,m)       !air water vapour defficit? Unit Pa
!                  co2ca=CO2treat*1.0E-6 ! CO2 concentration,ppm-->1.0E-6
!                  if(isite==2.and.yr>yrs_eq)then
!                      co2ca=(input_data(5,m)+200.)*1.0E-6
!                  endif

                  rain=input_data(5,m)    ! rain fal per hour
                  wind=ABS(input_data(6,m))     ! wind speed m s-1
                  PAR=input_data(7,m)             ! Unit ? umol/s/m-2
                  radsol=input_data(7,m)        ! unit ? PAR actually
!                  Rnet=input_data(9,m)
!!       endof Duke Forest

!                 Ajust some unreasonable values
                  RH=AMAX1(0.01,AMIN1(99.99,RH))
                  eairP = esat(Tair)*RH/100.             ! Added for SPRUCE, due to lack of VPD data
                  Dair=esat(Tair)-eairP
                  radsol=AMAX1(radsol,0.01)
                  hoy=hoy+1

                m=m+1											


                  if(radsol.gt.10.0) then
                      G=-25.0
                  else
                      G=20.5											
                  endif
                  Esoil=0.05*radsol
                  if(radsol.LE.10.0) Esoil=0.5*G						
                  Hcrop=0.1  ! never used in routine
                  Ecstot=0.1 ! never used in routine
                  Anet=0.1 ! never used in routine
                  DepH2O=0.2
!                 for daily mean conditions
                  VPD_d=VPD_d+Dair/24./1000.               
                  PAR_D=PAR_D+radsol/dtimes                 ! umol photons m-2 s-1   
                  ta= ta + tair/24.0             ! sum of a day, for calculating daily mean temperature
                  Ts=Ts+Tsoil/24.0
                  rain_d=rain_d+rain
!                 calculating scaling factor of NSC
                  if(NSC.le.NSCmin)fnsc=0.0
                  if(NSC.ge.NSCmax)fnsc=1.0
                  if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                     fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
                  endif
!                 update vcmx0 and eJmx0 according to C/N of leaves
                  Vcmx0 = Vcmax0*SNvcmax*1.0e-6
!                  eJmx0 = 2.7*Vcmx0  ! original
                  eJmx0 = 1.67*Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002
                  
                  call canopy(gpp,evap,transp,Acanop,Hcanop,Rsoilabs,  & ! outputs
     &              fwsoil,topfws, &                   ! from soil model
     &              LAI,Sps,&
     &              doy,hour,radsol,tair,dair,eairP,    &        ! from climate data file,including 
     &              wind,rain,&
     &              Rnet,G,Esoil,Hcrop,Ecstot,Anet,&
     &              Tsoil,DepH2O,&
     &              wsmax,wsmin,  &                              !constants specific to soil and plant
     &              lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,&
     &              stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &              Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
     &              H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,&
     &              conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
     &              Edjm,Entrpy,gam0,gam1,gam2,wcl,gddonset)

            call soilwater(wsmax,wsmin,rdepth,FRLEN,THKSL,    &   !constants specific to soil/plant
     &                rain,transp,evap,wcl,runoff,infilt,     &   !inputs
     &                fwsoil,topfws,omega)                      !outputs
                  ET=evap+transp
                  rain_yr=rain_yr+rain
                  transp_yr=transp_yr+transp
                  evap_yr=evap_yr+evap
                  runoff_yr=runoff_yr+runoff

                  call respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,&
     &                       Q10,Rl0,Rs0,Rr0,SNRauto,&
     &                       LAI,SLA,bmstem,bmroot,bmleaf,&
     &                       StemSap,RootSap,NSC,fnsc,&
     &                       RmLeaf,RmStem,RmRoot,Rmain)

!                 THE Third Part: update LAI
                  call plantgrowth(Tair,omega,GLmax,GRmax,GSmax,&
     &                    LAI,LAIMAX,LAIMIN,SLA,TauC(1),         &    !Tau_L,
     &                    bmleaf,bmroot,bmstem,bmplant,&
     &                    Rootmax,Stemmax,SapS,SapR,&
     &                    StemSap,RootSap,Storage,GDD5,&
     &                    stor_use,onset,accumulation,gddonset,&
     &                    Sps,NSC,fnsc,NSCmin,NSCmax,&
     &                    NSN,CN,CN0,SNgrowth,N_deficit,&
     &                    store,add,L_fall,ht,&
     &                    NPP,alpha_L,alpha_W,alpha_R,&
     &                    RgLeaf,RgStem,RgRoot,Rgrowth)

!                 THE Fourth PART: simulating C influx allocation in pools
                  call TCS_CN(Tair,Tsoil,omega,runoff,&
     &               NPP,alpha_L,alpha_W,alpha_R,L_fall,&
     &               tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
     &               CNmin,CNmax,NSNmax,NSNmin,alphaN,   &         ! nitrogen
     &               NSN,N_uptake,N_miner,QN,QNminer,&
     &               CN,CN0,fnsc,rdepth,N_deficit,&
     &               N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
     &               N_deposit,N_fixation,N_leach,N_vol,&
     &               SNvcmax,SNgrowth,SNRauto,SNrs,Q10)

!                 update NSC
                  Rauto  =Rmain+Rgrowth+Rnitrogen
                  NSC    =NSC+GPP-Rauto-(NPP-add)-store
                  Difference = GPP-Rauto-NPP
                  if(NSC<0)then
                      bmstem=bmstem+NSC/0.48
                      NPP=NPP+NSC
                      NSN=NSN-NSC/CN(2)
                      NSC=0.
                  endif
                  GL_d   =GL_d+NPP*alpha_L
                  GW_d   =GW_d+NPP*alpha_W
                  GR_d   =GR_d+NPP*alpha_R
                  LFALL_d=LFALL_d+L_fall
!                 update
                  RaLeaf = RgLeaf + RmLeaf
		  RaStem = RgStem + RmStem
		  RaRoot = RgRoot + RmRoot + Rnitrogen
                  WFALL_d=WFALL_d+OutC(2) !_wood
                  RFALL_d=RFALL_d+OutC(3) !_root
                  N_LG_d=N_LG_d+N_leaf
                  N_WG_d=N_WG_d+N_wood
                  N_RG_d=N_RG_d+N_root
                  N_LF_d=N_LF_d+N_LF
                  N_WF_d=N_WF_d+N_WF
                  N_RF_d=N_RF_d+N_RF

                  N_up_d=N_up_d+N_uptake
                  N_fix_d=N_fix_d+N_fixation
                  N_dep_d=N_dep_d+N_deposit
                  N_leach_d=N_leach_d+N_leach
                  N_vol_d=N_vol_d+N_vol

                  N_up_yr=N_up_yr+N_uptake
                  N_fix_yr=N_fix_yr+N_fixation
                  N_dep_yr=N_dep_yr+N_deposit
                  N_leach_yr=N_leach_yr+N_leach
                  N_vol_yr=N_vol_yr+N_vol

                  R_Ntr_yr=R_Ntr_yr + Rnitrogen
                  
           ! ==================== test variables
                  topfws_yr = topfws_yr+topfws/8760.
                  omega_yr=omega_yr+omega/8760.
                  fwsoil_yr=fwsoil_yr+fwsoil/8760.

!                 Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
                  Rhetero= Rh_pools(1)+Rh_pools(2)+Rh_pools(3) &
     &                    +Rh_pools(4)+Rh_pools(5)
                  Rsoil  =Rhetero+RmRoot+RgRoot+Rnitrogen
                  NEE=Rauto+Rhetero - GPP
                  Q_soil=QC(6) + QC(7) + QC(8)

                  bmleaf=QC(1)/0.48
                  bmstem=QC(2)/0.48
                  bmroot=QC(3)/0.48
                  bmplant=bmleaf+bmroot+bmstem
                  LAI=bmleaf*SLA
                  NMIN_d = NMIN_d+N_miner
!                 output hourly
                  Recoh=Rhetero+Rauto
                  ETh =ET !*1000.
                  Th  =transp !*1000.
                  Eh  =evap !*1000.
                  INTh=-9999
                  VPDh=Dair/1000.
                  ROh =runoff !*1000.
                  DRAINh=-9999
                  LEh =ETh*((2.501-0.00236*Tair)*1000.0)/3600.
                  SHh =-9999
                  LWh =-9999
                  NEP=-NEE

!                 sums of a day
                  diff_d=diff_d+difference
                  gpp_d=gpp_d + GPP
                  gpp_ra=gpp_ra+Rauto
                  NPP_d   =NPP_d+NPP
                  NEP_d=NEP_d+NEP
                  NEE_d=NEE_d+NEE
                  RECO_d=RECO_d+Recoh
                  Rh_d=  Rh_d + Rhetero
                  Ra_d=Reco_d-Rh_d
                  RLEAV_d=RLEAV_d+RmLeaf+RgLeaf
                  RWOOD_d=RWOOD_d+RmStem+RgStem
                  RROOT_d=RROOT_d+RmRoot+RgRoot+Rnitrogen
                  Rsoil_d=Rh_d+RROOT_d
                  NUP_d=NUP_d+N_uptake
                  NVOL_d=NVOL_d+N_vol
                  NLEACH_d=NLEACH_d+N_leach
                  transp_d=transp_d + transp*(24./dtimes)
                  evap_d=evap_d + evap*(24./dtimes)
                  ET_d=transp_d + evap_d
                  LE_d=LE_d+LEh/24.
                  Hcanop_d=Hcanop_d+Hcanop/(24./dtimes)
                  runoff_d=runoff_d+runoff
!                 sum of the whole year
                  diff_yr = diff_yr+difference
                  gpp_yr=gpp_yr+gpp
                  NPP_yr=NPP_yr+NPP
                  Rh_yr =Rh_yr +Rhetero
                  Ra_yr=Ra_yr+Rauto
                  Rh4_yr=Rh4_yr+Rh_pools(1)
                  Rh5_yr=Rh5_yr+Rh_pools(2)
                  Rh6_yr=Rh6_yr+Rh_pools(3)
                  Rh7_yr=Rh7_yr+Rh_pools(4)
                  Rh8_yr=Rh8_yr+Rh_pools(5)
                  Pool1 = Pool1+QC(1)/8760.
                  Pool2 = Pool2+QC(2)/8760.
                  Pool3 = Pool3+QC(3)/8760.
                  Pool4 = Pool4+QC(4)/8760.
                  Pool5 = Pool5+QC(5)/8760.
                  Pool6 = Pool6+QC(6)/8760.
                  Pool7 = Pool7+QC(7)/8760.
                  Pool8 = Pool8+QC(8)/8760.
                  out1_yr=out1_yr+OutC(1)
                  out2_yr=out2_yr+OutC(2)
                  out3_yr=out3_yr+OutC(3)
                  out4_yr=out4_yr+OutC(4)
                  out5_yr=out5_yr+OutC(5)
                  out6_yr=out6_yr+OutC(6)
                  out7_yr=out7_yr+OutC(7)
                  out8_yr=out8_yr+OutC(8)
                  NEE_yr=NEE_yr+NEE
                  GL_yr=GL_yr+NPP*alpha_L
                  GW_yr=GW_yr+NPP*alpha_W
                  GR_yr=GR_yr+NPP*alpha_R
!                 numbering         
                  n=n+1

                  if(isnan(gpp))then
                      write(*,*)'gpp is nan'
                      return
                  endif

              enddo              ! end of dtimes
              if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
                pheno=days
                phenoset=1
              endif

!             +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              
!             output results of canopy and soil models, daily
!             daily output
              INT_d=-9999
              DRAIN_d=-9999
              SH_d=-9999
              CSOIL_d=QC(6)+QC(7)+QC(8)
              LMA_d=QC(1)/0.48
              NSOIL_d=QN(6)+QN(7)+QN(8)+QNminer

              Rgrowth_d=-9999
              abvLitter=QC(4)  !-9999
              blLitter=QC(5) !-9999
!            write(*,*)yr,days,LAI,gpp_d,npp_d

!            if(yr.gt.yrs_eq)then  
                daily=daily+1
  
            Simu_dailyflux(1,daily)=GPP_d ! Leaf
            Simu_dailyflux(2,daily)=NEE_d	! Wood
            Simu_dailyflux(3,daily)=Reco_d	! Coarse roots
            Simu_dailyflux(4,daily)=QC(1)!*1.5
            Simu_dailyflux(5,daily)=GL_yr!*1.5
            Simu_dailyflux(6,daily)=QC(2)!*0.48
            Simu_dailyflux(7,daily)=GW_yr!*0.48
            Simu_dailyflux(8,daily)=QC(3)
            Simu_dailyflux(9,daily)=GR_yr
            Simu_dailyflux(10,daily)=(QC(6)+QC(7)+QC(8))!*13.8   ! Soil
            Simu_dailyflux(11,daily)=pheno   ! Soil
            Simu_dailyflux(12,daily)=LAI !QC(1)/(QC(1)+QC(4)) 
            
!            endif
!            if(yr.ge.(yrlim-first_year+1) .and. days.ge.dylim) goto 650
        enddo                         ! end of idays
                
            storage=accumulation
            stor_use=Storage/times_storage_use
            if(yr.eq.yrs_eq+yr_length .and. MCMC.eq.1)then
            write(*,*)yr,LAI,gpp_yr,NPP_yr,pheno
            write(61,601)year,LAI,gpp_yr,NPP_yr,real(pheno)
            endif
                
            if(MCMC.ne.1) then
                write(*,*)year,LAI,gpp_yr,NPP_yr,pheno
                write(61,601)year,LAI,gpp_yr,NPP_yr,Ra_yr,Rh_yr, &
                &   ET,rain_yr,transp_yr,evap_yr,runoff_yr,GL_yr,    &
                &   GW_yr,GR_yr,Pool1,Pool2,Pool3,Pool4,Pool5,   &
                &   Pool6,Pool7,Pool8,out1_yr,out2_yr,out3_yr,   &
                &   out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
            endif            
601         format(i7,",",29(f15.4,","))
            accumulation=0.0
            onset=0
         enddo            !end of simulations multiple years
         
650      if(MCMC.ne.1)then
         do i=1,daily
         write(62,602)i,(Simu_dailyflux(j,i),j=1,12)
         enddo
       
         endif
602     format((i7),",",12(f15.4,","))
603     format((i7),",",10(f15.4,","))
999      continue

    return
    end
      

!     ****************************************************************************
      subroutine consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,                &
     &   Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,    &
     &   wleaf,gsw0,Vcmx0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,       &
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2) 
     
      real tauL(3), rhoL(3), rhoS(3)
      pi = 3.1415926
!     physical constants
      tauL(1)=0.1                  ! leaf transmittance for vis
      rhoL(1)=0.1                  ! leaf reflectance for vis
      rhoS(1)=0.1                  ! soil reflectance for vis
      tauL(2)=0.425                ! for NIR
      rhoL(2)=0.425                ! for NIR
      rhoS(2)=0.3                  ! for NIR - later function of soil water content
      tauL(3)=0.00                 ! for thermal
      rhoL(3)=0.00                 ! for thermal
      rhoS(3)=0.00                 ! for thermal
      emleaf=0.96
      emsoil=0.94
      Rconst=8.314                 ! universal gas constant (J/mol)
      sigma=5.67e-8                ! Steffan Boltzman constant (W/m2/K4)
      cpair=1010.                  ! heat capapcity of air (J/kg/K)
      Patm=101325. !1.e5           ! atmospheric pressure  (Pa)
      Trefk=293.2                  !reference temp K for Kc, Ko, Rd
      H2OLv0=2.501e6               !latent heat H2O (J/kg)
      AirMa=29.e-3                 !mol mass air (kg/mol)
      H2OMw=18.e-3                 !mol mass H2O (kg/mol)
      chi=0.93                     !gbH/gbw
      Dheat=21.5e-6                !molecular diffusivity for heat
!     plant parameters
      gsw0 = 1.0e-2                !g0 for H2O in BWB model
      eJmx0 = Vcmx0*2.7            !@20C Leuning 1996 from Wullschleger (1993)
      theta = 0.9
      wleaf=0.01                   !leaf width (m)

!     thermodynamic parameters for Kc and Ko (Leuning 1990)
      conKc0 = 302.e-6                !mol mol^-1
      conKo0 = 256.e-3                !mol mol^-1
      Ekc = 59430.                    !J mol^-1
      Eko = 36000.                    !J mol^-1
!     Erd = 53000.                    !J mol^-1
      o2ci= 210.e-3                   !mol mol^-1

!     thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
      Eavm = 116300.               !J/mol  (activation energy)
      Edvm = 202900.               !J/mol  (deactivation energy)
      Eajm = 79500.                !J/mol  (activation energy) 
      Edjm = 201000.               !J/mol  (deactivation energy)
      Entrpy = 650.                !J/mol/K (entropy term, for Jmax & Vcmax)

!     parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
      gam0 = 28.0e-6               !mol mol^-1 @ 20C = 36.9 @ 25C
      gam1 = .0509
      gam2 = .0010
      return
      end
!****************************************************************************

!      a sub-model for calculating C flux and H2O flux of a canopy
!      adapted from a two-leaf canopy model developed by Wang Yingping
         
        subroutine canopy(gpp,evap,transp,Acanop,Hcanop,Rsoilabs, &  ! outputs
     &                    fwsoil,topfws,           &! from soil model
     &                    LAI,Sps,&
     &               doy,hour,radsol,tair,Dair,eairP,&! from climate data file,including 
     &               wind,rain,&
     &               Rnet,G,Esoil,Hcrop,Ecstot,Anet,&
     &               Tsoil,DepH2O,&
     &               wsmax,wsmin,&  !constants specific to soil and plant
     &               lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,&
     &               stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
     &               Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
     &               H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,&
     &               conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
     &               Edjm,Entrpy,gam0,gam1,gam2,wcl,gddonset)

        real lat,doy
      real gpp,evap,transp,LAI,Rsoilabs
      real tauL(3),rhoL(3),rhoS(3),reffbm(3),reffdf(3)
      real extkbm(3),extkdm(3)
      real Radabv(2),Qcan(3,2)
      real gddonset
!     extra variables used to run the model for the wagga data
      real topfws        ! from siol subroutine      
      integer idoy,ihour,ileaf
      integer jrain,i,j,k

!     additional arrays to allow output of info for each layer
      real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      real GbwcL(5),GswcL(5),hG(5),hIL(5)
      real Gaussx(5),Gaussw(5),Gaussw_cum(5)
      real wcl(10)
      character*80 commts
!     Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!     5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
      data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/

!     calculate beam fraction in incoming solar radiation
      call  yrday(doy,hour,lat,radsol,fbeam)
      idoy=int(doy)
      hours=idoy*1.0+hour/24.0
      coszen=sinbet(doy,lat,pi,hour)             !cos zenith angle of sun

!     set windspeed to the minimum speed to avoid zero Gb
      if(wind.lt.0.01) wind=0.01
!     calculate soil albedo for NIR as a function of soil water (Garratt pp292)
      if(topfws.gt.0.5) then
            rhoS(2)=0.18
      else
            rhoS(2)=0.52-0.68*topfws
      endif
!        assign plant biomass and leaf area index at time t
!        assume leaf biomass = root biomass
      FLAIT =LAI 
      eairP=esat(Tair)-Dair                !air water vapour pressure
      radabv(1)=0.5*radsol                 !(1) - solar radn
      radabv(2)=0.5*radsol                 !(2) - NIR
!     call multilayer model of Leuning - uses Gaussian integration but radiation scheme
!     is that of Goudriaan
      call xlayers(Sps,Tair,Dair,radabv,G,Esoil,fbeam,eairP,&
     &           wind,co2ca,fwsoil,wcl,LAI,coszen,idoy,hours,&
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,&
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,&
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,&
     &           gsw0,alpha,stom_n,wsmax,wsmin,&
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,&
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,&
     &           extKb,Rsoilabs,Hsoil,Acan1,Acan2,Ecan1,Ecan2,&
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset)


         Acanop=Acan1+Acan2
         Ecanop=Ecan1+Ecan2
         gpp=Acanop*3600.0*12.0                           ! every hour, g C m-2 h-1
         transp=AMAX1(Ecanop*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.) ! mm H2O /hour
         evap=AMAX1(Esoil*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.)
!        H2OLv0=2.501e6               !latent heat H2O (J/kg)

      return
      end
      
!****************************************************************************
!   autotrophic respiration
    subroutine respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,&
     &                       Q10,Rl0,Rs0,Rr0,SNRauto,&
     &                       LAI,SLA,bmstem,bmroot,bmleaf,&
     &                       StemSap,RootSap,NSC,fnsc,&
     &                       RmLeaf,RmStem,RmRoot,Rmain)
!     calculate plant and soil respiration by the following equation:
!     RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
    implicit none
    real LAIMIN,LAI,GPP,SLA
    real Tair,Tsoil,DepH2O
    real bmstem,bmroot,bmleaf,StemSap,RootSap
    real NSC,fnsc
    real Q10
    real RmLeaf,RmStem,RmRoot,Rmain
    real Rl0,Rs0,Rr0,SNRauto
    real conv                  ! converter from "umol C /m2/s" to "gC/m2/hour"

    conv=3600.*12./1000000.    ! umol C /m2/s--> gC/m2/hour
    if(LAI.gt.LAIMIN) then
        RmLeaf=Rl0*SNRauto*bmleaf*0.48*SLA*0.1     &               
            &   *Q10**((Tair-10.)/10.)*fnsc*conv
        RmStem=Rs0*SNRauto*StemSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
        RmRoot=Rr0*SNRauto*RootSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
    else
        RmLeaf=0.3*GPP
        RmStem=0.3*GPP
        RmRoot=0.4*GPP
    endif
        Rmain=Rmleaf+Rmstem+Rmroot
    if(Rmain > 0.0015*NSC)then             ! If Total autotropic respiration greater than 0.15% of Nonstructure Carbon, rescale. 
        Rmleaf=Rmleaf/Rmain*0.0015*NSC
        Rmstem=Rmstem/Rmain*0.0015*NSC
        Rmroot=Rmstem/Rmain*0.0015*NSC
        Rmain=Rmleaf+Rmstem+Rmroot
    endif
    return
    end

!*******************************************************************
!     subroutine for soil moisture
    subroutine soilwater(wsmax,wsmin,rdepth,FRLEN,THKSL,    &   !constants specific to soil/plant
     &                rain,transp,evap,wcl,runoff,infilt,   &   !inputs
     &                fwsoil,topfws,omega)                      !outputs
!     All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
    implicit none
!   soil traits
    real wsmax,wsmin,wsmaxL(10),wsminL(10) !from input percent x%
    real FLDCAP,WILTPT ! ie. 0.xx
!   plant traits
    real rdepth
    integer nfr
!   climate conditions
    real rain ! mm/hour
!   output from canopy model
    real evap,transp
!   output variables
    real fwsoil,topfws,omega
    real fw(10),ome(10)

    real thksl(10),depth(10),wsc(10),WUPL(10),EVAPL(10),SRDT(10)
    real plantup(10)
    real Tsrdt
    real frlen(10) !fraction of root length in every layer
    real wcl(10) !volum ratio
 !      real fwcln(10) !  fraction of water in layers, like field capacity
    real DWCL(10),Tr_ratio(10)
    real wtadd,twtadd,infilt,runoff,tr_allo

    real exchangeL,supply,demand,omegaL(10)
    integer i,j,k
    real infilt_max

    infilt_max=15.
    WILTPT =wsmin/100.0
    FLDCAP =wsmax/100.0
    
    do i=1,10
        dwcl(i)=0.0
        evapl(i)=0.0
        WUPL(i)=0.0
        SRDT(i)=0.0
        DEPTH(i)=0.0
    enddo

!   Determine which layers are reached by the root system. 
!   Layer volume (cm3)
    DEPTH(1)=10.0
    DO i=2,10
        DEPTH(i)=DEPTH(i-1)+THKSL(i)
    enddo
    do i=1,10
        IF(rdepth.GT.DEPTH(i)) nfr=i+1
    enddo
    IF (nfr.GT.10) nfr=10

! *** water infiltration through layers
    infilt=infilt+rain  !mm/hour
!   Loop over all soil layers.
    TWTADD=0
    IF(infilt.GT.0.0)THEN
!       Add water to this layer, pass extra water to the next.
        WTADD=AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
!       change water content of this layer
        WCL(1)=(WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
!       FWCLN(I)=WCL(I)       !  /VOLUM(I)! update fwcln of this layer
        TWTADD=TWTADD+WTADD       !calculating total added water to soil layers (mm)
        INFILT=INFILT-WTADD !update infilt
    ENDIF

!   runoff
    runoff=INFILT*0.0019   ! Shuang Modifed
    infilt = infilt-runoff
    
!   water redistribution among soil layers
    do i=1,10
        wsc(i)=Amax1(0.00,(wcl(i)-wiltpt)*THKSL(i)*10.0)
        omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAP-WILTPT))
    enddo
    supply=0.0
    demand=0.0
    do i=1,9
        if(omegaL(i).gt.0.3)then
            supply=wsc(i)*omegaL(i)
            demand=(FLDCAP-wcl(i+1))*THKSL(i+1)*10.0      &
                &               *(1.0-omegaL(i+1))
            exchangeL=AMIN1(supply,demand)
            wsc(i)=wsc(i)- exchangeL
            wsc(i+1)=wsc(i+1)+ exchangeL
            wcl(i)=wsc(i)/(THKSL(i)*10.0)+wiltpt
            wcl(i+1)=wsc(i+1)/(THKSL(i+1)*10.0)+wiltpt
        endif
    enddo
    wsc(10)=wsc(10)-wsc(10)*0.00001     ! Shuang modifed
    runoff = runoff+wsc(10)*0.00001     ! Shuang modifed
    wcl(10)=wsc(10)/(THKSL(10)*10.0)+wiltpt
!    end of water redistribution among soil layers

!   Redistribute evaporation among soil layers
    Tsrdt=0.0
    DO i=1,10
!   Fraction of SEVAP supplied by each soil layer
    SRDT(I)=EXP(-6.73*(DEPTH(I)-THKSL(I)/2.0)/100.0) !/1.987
!   SRDT(I)=AMAX1(0.0,SRDT(I)*(wcl(i)-wiltpt)) !*THKSL(I))
    Tsrdt=Tsrdt+SRDT(i)  ! to normalize SRDT(i)
    enddo
    do i=1,10
        EVAPL(I)=Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)  !mm
        DWCL(I)=EVAPL(I)/(THKSL(I)*10.0) !ratio
        wcl(i)=wcl(i)-DWCL(i)
    enddo
    evap=0.0       
    do i=1,10
        evap=evap+EVAPL(I)
    enddo

!   Redistribute transpiration according to root biomass
!   and available water in each layer
    tr_allo=0.0
    do i=1,nfr
        tr_ratio(i)=FRLEN(i)*wsc(i) !*(wcl(i)-wiltpt)) !*THKSL(I))
        tr_allo=tr_allo+tr_ratio(i)
    enddo
    do i=1,nfr
        plantup(i)=AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm              
        wupl(i)=plantup(i)/(thksl(i)*10.0)
        wcl(i)=wcl(i)-wupl(i)
    enddo
    transp=0.0
    do i=1,nfr
        transp=transp+plantup(i)
    enddo
    
!   Output fwsoil, omega, and topfws
    do i=1,nfr       
        ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
        ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
        fw(i)=amin1(1.0,3.333*ome(i))
    enddo
    topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))
    fwsoil=0.0
    omega=0.0
    do i=1,nfr
        fwsoil=fwsoil+fw(i)*frlen(i)
        omega=omega+ome(i)*frlen(i)
    enddo
    return
    end
    
!**********************************************************************
!     plant growth model
    subroutine plantgrowth(Tair,omega,GLmax,GRmax,GSmax,&
     &                       LAI,LAIMAX,LAIMIN,SLA,Tau_L,&
     &                       bmleaf,bmroot,bmstem,bmplant,&
     &                       Rootmax,Stemmax,SapS,SapR,&
     &                       StemSap,RootSap,Storage,GDD5,&
     &                       stor_use,onset,accumulation,gddonset,&
     &                       Sps,NSC,fnsc,NSCmin,NSCmax,&
     &                       NSN,CN,CN0,SNgrowth,N_deficit,&
     &                       store,add,L_fall,ht,&
     &                       NPP,alpha_L,alpha_W,alpha_R,&
     &                       RgLeaf,RgStem,RgRoot,Rgrowth)
      implicit none
      real NSC,NSCmin,NSCmax,fnsc,N_deficit
      real CN(8),CN0(8),NSN,nsCN
      real SnscnL,SnscnS,SnscnR
      real store,Storage,GDD5,stor_use,accumulation,gddonset
      integer onset
      real GLmax,GRmax,GSmax,TauLeaf
      real GrowthP,GrowthL,GrowthR,GrowthS
      real Tair,omega,LAI,LAIMAX,LAIMIN,SLA
!     biomass
      real bmleaf,bmroot,bmstem,bmplant,NPP
      real ht,hmax,hl0,CNP0
      REAL LAIMAX0,la0,GPmax,acP,c1,c2
      real Rootmax,Stemmax,SapS,SapR
      real bmL,bmR,bmP,bmS,StemSap,RootSap
      real Rgrowth,Rgroot,Rgleaf,Rgstem
!     scalars
      real St,Sw,Ss,Sn,SL_rs,SR_rs,Slai,Sps,SNgrowth,phiN
      real RS,RS0,RSw
      real gamma_W,gamma_Wmax,gamma_T,gamma_Tmax,gamma_N
      real beta_T,Tcold,Twarm,Topt
      real bW,bT,W
      real L_fall,L_add,add,NL_fall,NL_add,Tau_L
      real alpha_L,alpha_W,alpha_R,alpha_St
      integer i

    Twarm=35.0
    Tcold=5.0
!    Tcold=0.0       ! For SPRUCE
    Topt=30.
    phiN=0.33

    bmL=bmleaf*0.48   ! Carbon
    bmR=bmRoot*0.48
    bmS=bmStem*0.48

    if(bmL.lt.NSC/0.333)bmL=NSC/0.333
    if(bmR.lt.NSC/0.333)bmR=NSC/0.333
    if(bmS.lt.NSC/0.334)bmS=NSC/0.334
    StemSap=SapS*bmS  ! Weng 12/05/2008
    RootSap=SapR*bmR
    if(StemSap.lt.0.001)StemSap=0.001
    if(RootSap.lt.0.001)RootSap=0.001

    bmP=bmL+bmR+bmS					! Plant C biomass 
    acP=bmL+StemSap+bmS					! Plant available sapwood C  
    CNp0=bmP/(bmL/CN0(1)+bmR/CN0(3)+bmS/CN0(2))		! Plant CN ratio

    hmax=24.19   ! m
    hl0=0.00019  ! m2/kg C
    LAIMAX0=6.
    la0=0.2
    ht=hmax*(1.-exp(-hl0*bmP))				! Scaling plant C biomass to height
    LAIMAX=AMAX1(LAIMAX0*(1.-exp(-la0*ht)),LAIMIN+0.1)  ! Scaling plant height to maximum LAI

!   Phenology
    if((GDD5.gt.gddonset).and.onset.eq.0.and.storage.gt.stor_use) then
        onset=1
    endif
    if((onset.eq.1).and.(storage.gt.stor_use))then
        if(LAI.lt.LAIMAX)add=stor_use
        storage=storage-add
    else
        add=0.0
        onset=0
    endif
    if(accumulation.lt.(NSCmax+0.005*RootSap))then
        store=AMAX1(0.,0.005*NSC)			! 0.5% of nonstructure carbon is stored
    else
        store=0.0
    endif
    accumulation=accumulation+store

!   Scalars for plant growth
!      Sps=Amin1(1.0,3.33*AMAX1(0.0,1.0 - fnsc))
      Sps=Sps*(1.-exp(-phiN*NSN))							! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
      Ss=AMIN1(1.0,2.*fnsc)
      RS0=1.0
      RS=bmR/bmL
      SL_rs=RS/(RS+RS0*(2.-W))
      SR_rs=(RS0*(2.-W))/(RS+RS0*(2.-W))
      Slai=amin1(1.0,2.333*(LAIMAX-LAI)/(LAIMAX-LAIMIN))
      St=AMAX1(0.0, 1.0-exp(-(Tair-gddonset/10.)/5.0))  !0.5 !
!      Sw=AMAX1(0.333, 0.333+omega)
      Sw=AMIN1(0.5, AMAX1(0.333, 0.333+omega))
      W = AMIN1(1.0,3.333*omega)

!     Plant growth and allocation, based on LM3V
      GPmax=(GLmax*bmL+GSmax*StemSap+GRmax*bmR) !/acP					
      GrowthP=AMIN1(GPmax*fnsc*St*(1.-exp(-NSN)),  & ! 
     &              0.004*NSC,&
     &              0.004*NSN*CNp0)
 
     !      c1=(bmR+200.)/bmL*CN(1)/CN0(1) !+N_deficit/NSN
!      c1=bmL/bmR*CN(1)/CN0(1) !+N_deficit/NSN
!      c2=0.5*250e3*SLA*0.00021*ht*2.
      !write(*,*)LAI
      !GrowthL=MAX(0.0,MIN(GrowthP*0.5,0.05*(LAIMAX-LAI)/SLA))    ! 1./(1.+c1+c2)
      GrowthL=MAX(0.0,GrowthP*0.5) 
      GrowthR=MIN(GrowthP*0.4,MAX(0.0,0.75/Sw*bmL-bmR))  ! *c1/(1.+c1+c2)
      GrowthS=MAX(0.0,GrowthP - (GrowthL+GrowthR) )         ! *c2/(1.+c1+c2)

      NPP = GrowthL + GrowthR + GrowthS + add       ! Modified by Jiang Jiang 2015/10/13

      if(NPP.eq.0.0)then
            alpha_L=0.333
            alpha_W=0.333
            alpha_R=0.333
      else
            alpha_L=(GrowthL+add)/NPP     
            alpha_W=GrowthS/NPP
            alpha_R=GrowthR/NPP
      endif
!     Carbon cost for growth
!     Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
      Rgleaf=0.5*GrowthL
      Rgstem=0.5*GrowthS
      Rgroot=0.5*GrowthR
      Rgrowth=Rgleaf+Rgstem+Rgroot

!     Leaf litter 

      gamma_Wmax=0.12/24. ! maxmum leaf fall rate per hour
      gamma_Tmax=0.12/24.

      bW=4.0
      bT=2.0

      if(Tair.gt.(Tcold+10.)) then
            beta_T=1.
      else 
            if(Tair.gt.Tcold)beta_T=(Tair-Tcold)/10.
            if(Tair.LE.Tcold)beta_T=0.0
      endif

      if (tau_L < 8760.)then
                gamma_W=(1. - W)     **bW  * gamma_Wmax
                gamma_T=(1. - beta_T)**bT * gamma_Tmax
      else
                gamma_W=0.
                gamma_T=0.
      endif
      gamma_N=1.0/Tau_L*Sw      ! Modify by Jiang Jiang 2015/10/20
      if(LAI < LAIMIN) then
            gamma_W=0.
            gamma_T=0.
            gamma_N=0.
      endif
    !  L_fall=bmleaf*0.48*AMIN1((gamma_T+gamma_N),0.99)
      L_fall=bmleaf*0.48*gamma_N

      return
      end

!************************************************************************
!     carbon transfer according to Xu et al. 2007 
      subroutine TCS_CN(Tair,Tsoil,omega,runoff,&
     &               NPP,alpha_L,alpha_W,alpha_R,L_fall,&
     &               tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
     &               CNmin,CNmax,NSNmax,NSNmin,alphaN,    &        ! nitrogen
     &               NSN,N_uptake,N_miner,QN,QNminer,&
     &               CN,CN0,fnsc,rdepth,N_deficit,&
     &               N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
     &               N_deposit,N_fixation,N_leach,N_vol,&
     &               SNvcmax,SNgrowth,SNRauto,SNrs,Q10) 
         
      implicit none
      real NPP,NPP_L,NPP_W,NPP_R
      real L_fall,L_add,LAI,SLA,rdepth
      real Tair,Tsoil,omega,runoff
!     allocation ratios
      real alpha_L,alpha_W,alpha_R
!     pools
      real Q_plant,QC(8),TauC(8),OutC(8)
      real etaL,etaW,etaR                ! the percentage of fine litter of the litters from plant parts
      real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M
      real Rh_pools(5),Q10h(5),Q10 ! Q10 of the litter and soil C pools
!     the fraction of C-flux which enters the atmosphere from the kth pool
      real f_CO2_fine,f_CO2_coarse,f_CO2_Micr,f_CO2_Slow,f_CO2_Pass
!     for nitrogen sub-model
      real CN0(8),CN(8),OutN(8),QN(8),QNminer,QNplant
      real CNmin,CNmax,NSNmax,NSNmin,NSN
      real CN_plant,CN_foliage
      real N_demand,N_deficit,N_immob,N_imm(5),N_fixation,Nfix0
      real N_transfer,N_miner,N_uptake,N_deposit,N_loss,N_leach,N_vol
      real alphaN,Qroot0,Cfix0
      real Scalar_N_flow,Scalar_N_T
      real N_leaf,N_wood,N_root
      real N_LF,N_WF,N_RF
      real NSC,fnsc,ksye
      real SNvcmax,SNgrowth,SNRauto,SNrs
      real kappaVcmax
      real SNfine,SNcoarse,SNmicr,SNslow,SNpass
      real Rnitrogen,costCuptake,costCfix,costCreuse
      real Creuse0,Nup0,N_deN0,LDON0
!     the variables relative to soil moisture calcualtion
      real S_omega    !  average values of the moisture scaling functions
      real S_t(5)     !  average values of temperature scaling functions
      real S_w_min    !  minimum decomposition rate at zero available water
!     For test
      real totalC1,totalN1,C_in,C_out,N_in,N_out,totalC2,totalN2
      real ScNloss

      integer i,j,k,n,m
      integer day,week,month,year

!     temperature sensitivity of Rh_pools
!      Q10h=(/2.0,2.0,2.0,2.0,2.0/)  ! for Oak Ridge
      Q10h=(/Q10,Q10,Q10,Q10,Q10/)
!
      Qroot0=500.
      Nfix0=1./60.   ! maximum N fix ratio, N/C
      Nup0 =0.02     ! nitrogen uptake rate
      Cfix0=12.      ! C cost per N for fixation
      ksye=0.05      ! C cost per N for uptake
      Creuse0=2.     ! C cost per N for resorption
      ScNloss=1.
      N_deN0=1.E-3*ScNloss   ! 1.E-3, 5.E-3, 10.E-3, 20.E-3
      LDON0=1.E-3*ScNloss

      Rnitrogen=0.
!     for N scalars
      CNmin=40.0
      CNmax=200.0
!     Max and min NSN pool
      NSNmax = QN(1) + 0.2*QN(2) + QN(3)  ! 15.0
      NSNmin = 0.01
!     partitioning coefficients
      etaL=0.6          ! 60% of foliage litter is fine, didn't use
      etaW=0.15        ! 15% of woody litter is fine
      etaR=0.85         ! 85% of root litter is fine  , didn't use    
      f_F2M=0.55        ! *exp((CN0(4)-CN_fine)*0.1)
      f_C2M=0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
      f_C2S=0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
      f_M2S=0.3
      f_M2P=0.1
      f_S2P=0.2        !0.03 Change by Jiang Jiang 10/10/2015
      f_S2M=0.5
      f_P2M=0.45

!     calculating soil scaling factors, S_omega and S_tmperature
      S_w_min=0.08 !minimum decomposition rate at zero soil moisture
      S_omega=S_w_min + (1.-S_w_min) * Amin1(1.0, 0.3*omega)

      do i=1,5
!        S_t(i)=Q10h(i)**((Tsoil-5.)/10.)  ! Oak
         S_t(i)=Q10h(i)**((Tsoil-10.)/10.)  ! Duke
      enddo

!     Calculating NPP allocation and changes of each C pool
      NPP_L = alpha_L * NPP           ! NPP allocation
      NPP_W = alpha_W * NPP
      NPP_R = alpha_R * NPP
!     N scalars on decomposition
      SNfine  =exp(-(CN0(4)-CN(4))/CN0(4)) 
      SNcoarse=exp(-(CN0(5)-CN(5))/CN0(5)) 
      SNmicr  =exp(-(CN0(6)-CN(6))/CN0(6)) 
      SNslow  =1. !exp(-(CN0(7)-CNC(7))/CN0(7)) 
      SNpass  =exp(-(CN0(8)-CN(8))/CN0(8)) 
      
!     the carbon leaving the pools
      OutC(1)=L_fall
      OutC(2)=QC(2)/tauC(2)*S_omega !*exp(CN(2)/CN0(2)-1.) 
      OutC(3)=QC(3)/tauC(3)*S_omega

      OutC(4)=QC(4)/tauC(4)*S_omega* S_T(1)*CN(4)/CN0(4)!*SNfine
      OutC(5)=QC(5)/tauC(5)*S_omega* S_T(2)*CN(5)/CN0(5)!*SNcoarse
      OutC(6)=QC(6)/tauC(6)*S_omega* S_T(3)!*SNmicr
      OutC(7)=QC(7)/tauC(7)*S_omega* S_T(4)!*SNslow
      OutC(8)=QC(8)/tauC(8)*S_omega* S_T(5)!*SNpass

!     heterotrophic respiration from each pool
      Rh_pools(1)=OutC(4)* (1. - f_F2M)
      Rh_pools(2)=OutC(5)* (1. - f_C2M - f_C2S)
      Rh_pools(3)=OutC(6)* (1. - f_M2S - f_M2P)
      Rh_pools(4)=OutC(7)* (1. - f_S2P - f_S2M)
      Rh_pools(5)=OutC(8)* (1. - f_P2M)

!========================================================================
!     Nitrogen part
!     nitrogen leaving the pools and resorption
      do i=1,8
         OutN(i) = OutC(i)/CN(i)
      enddo

!     nitrogen mineralization
      N_miner=OutN(4)* (1. - f_F2M)  &
     &       +OutN(5)* (1. - f_C2M - f_C2S) &
     &       +OutN(6)* (1. - f_M2S - f_M2P) &
     &       +OutN(7)* (1. - f_S2P - f_S2M) &
     &       +OutN(8)* (1. - f_P2M)

!     Nitrogen immobilization
      N_imm=0.
      N_immob=0.
      if(QNminer>0)then
          do i=4,8
            if(CN(i)>CN0(i))then
                N_imm(i-3)=Amin1(QC(i)/CN0(i)-QC(i)/CN(i)  &
     &             ,0.1*QNminer)
                N_immob=N_immob+N_imm(i-3)
            endif
          enddo
      endif

!     Let plant itself choose the strategy between using C to uptake
!     or fix N2 by comparing C invest.
!     N demand
      N_demand=NPP_L/CN0(1)+NPP_W/CN0(2)+NPP_R/CN0(3) !+N_deficit
!     Nitrogen input:
      N_transfer=0.
	N_uptake=0.
	N_fixation=0.
      costCuptake=0.
      costCfix=0.
      costCreuse=0.
!     1. Nitrogen resorption
      N_transfer=(OutN(1) + OutN(2) +OutN(3))*alphaN
      costCreuse= Creuse0*N_transfer
      N_demand=N_demand-N_transfer

      If(N_demand>0.0)then

!     2.  N uptake
          if(ksye/QNminer<Cfix0)then
              N_uptake=AMIN1(N_demand+N_deficit,      &
     &                       QNminer*QC(3)/(QC(3)+Qroot0), &
     &                       Nup0*NSC/(ksye/QNminer))
              costCuptake=N_uptake*(ksye/QNminer)
              N_demand=N_demand-N_uptake
          elseif(NSN<24.*30.*N_demand)then
!     3.  Nitrogen fixation
              N_fixation=Amin1(N_demand,fnsc*Nfix0*NSC)
              costCfix=Cfix0*N_fixation
              N_demand=N_demand-N_fixation
          endif
      endif
      N_deficit=N_deficit+N_demand

!     update NSN
      NSN=NSN+N_transfer+N_uptake+N_fixation
!     Total C cost for nitrogen
      Rnitrogen=costCuptake+costCfix+costCreuse

!      Oak Ridge N fixation rate: 
!      asymbiotic: 2 mg N/m2/yr ;  symbiotic: 65 mg/m2/yr, Oak Ridge
!      N_fixation=0.067/8760. ! Oak Ridge
!      N_fixation=0.23/8760.  ! Duke

!     Nitrogen using, non-structural nitrogen pool, NSN
      N_leaf =AMIN1(NPP*alpha_L/CN(1)+QC(1)/CN0(1)-QC(1)/CN(1),0.2*NSN)
      N_wood =AMIN1(NPP*alpha_W/CN(2)                         ,0.1*NSN)
      N_root =AMIN1(NPP*alpha_R/CN(3)+QC(3)/CN0(3)-QC(3)/CN(3),0.2*NSN)
      NSN=NSN-(N_leaf+N_wood+N_root)

      N_LF=OutN(1)*(1.-alphaN)
      N_WF=OutN(2)*(1.-alphaN)
      N_RF=OutN(3)*(1.-alphaN)

!     update QNminer
      QNminer=QNminer+N_miner+N_deposit  &
     &              -(N_uptake+N_immob)

!     Loss of mineralized N and dissolved organic N
      Scalar_N_flow=0.5*runoff/rdepth
!      Scalar_N_T=0.005*(Tsoil+273.)/(Tsoil+273+333.)
      Scalar_N_T=N_deN0*exp((Tsoil-25.)/10.)
      N_leach=Scalar_N_flow*QNminer+Scalar_N_flow*QN(6)*LDON0
      N_vol  =Scalar_N_T*QNminer
      N_loss =N_leach + N_vol

!     update QNminer
      QNminer=QNminer-N_loss


!     update plant carbon pools, ! daily change of each pool size
      QC(1)=QC(1) - OutC(1) + NPP_L    
      QC(2)=QC(2) - OutC(2) + NPP_W
      QC(3)=QC(3) - OutC(3) + NPP_R
      QC(4)=QC(4) - OutC(4) + OutC(1)+etaW*OutC(2)+OutC(3)
      QC(5)=QC(5) - OutC(5) + (1.-etaW)*OutC(2)
      QC(6)=QC(6) - OutC(6) + f_F2M*OutC(4)+f_C2M*OutC(5)     &         
     &                      + f_S2M*OutC(7)+f_P2M * OutC(8)
      QC(7)=QC(7) - OutC(7)+f_C2S*OutC(5)+f_M2S*OutC(6)
      QC(8)=QC(8) - OutC(8)+f_M2P*OutC(6)+f_S2P*OutC(7)

      Q_plant =QC(1) + QC(2) + QC(3)
!     update nitrogen pools
      QN(1)=QN(1) - OutN(1) + N_leaf
      QN(2)=QN(2) - OutN(2) + N_wood
      QN(3)=QN(3) - OutN(3) + N_root
      QN(4)=QN(4) - OutN(4)+ N_imm(1)     &
     &            + (OutN(1) + etaW*OutN(2) + OutN(3))*(1.-alphaN)
      QN(5)=QN(5) - OutN(5) + N_imm(2)   &
     &            + (1.-etaW)*OutN(2)*(1.-alphaN)

      QN(6)=QN(6) - OutN(6) + N_imm(3) - Scalar_N_flow*QN(6)*LDON0  &
     &            + f_F2M*OutN(4)+f_C2M*OutN(5)  &
     &            + f_S2M*OutN(7)+f_P2M*OutN(8)
      QN(7)= QN(7) - OutN(7) + N_imm(4)  &
     &                         + f_C2S*OutN(5) &
     &                         + f_M2S*OutN(6)
      QN(8)= QN(8) - OutN(8) + N_imm(5) &
     &         + f_M2P*OutN(6) + f_S2P*OutN(7)
      QNplant = QN(1) + QN(2)+ QN(3)

!     update C/N ratio
      CN=QC/QN
      CN_foliage=(QC(1)+QC(3))/(QN(1)+QN(3))

!     calculate N related scalars for Oak Ridge
!      SNvcmax =exp(-(CN(1)-CN0(1))) ! /CN0(1) ! Oak Ridge
!      SNgrowth=exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
!      SNRauto =AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
!      SNrs=1.

!     calculate N related scalars for Duke FACE
      kappaVcmax=CN0(1)/1.
      SNvcmax =exp(-kappaVcmax*(CN(1)-CN0(1))/CN0(1)) ! /CN0(1) ! Duke
      SNgrowth=exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
      SNRauto =exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
      SNrs=1.

      return
      end


!     ========================================================================================
!     subroutines used by canopy submodel
      subroutine xlayers(Sps,Tair,Dair,radabv,G,Esoil,fbeam,eairP,  &
     &           wind,co2ca,fwsoil,wcl,FLAIT,coszen,idoy,hours,   &
     &           tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,            &
     &           Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,           &
     &           cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,         &
     &           gsw0,alpha,stom_n,wsmax,wsmin,                     &
     &           Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &           Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &           extKb,Rsoilabs,Hsoil,Acan1,Acan2,Ecan1,Ecan2,      &
     &           RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset)


!    the multi-layered canopy model developed by 
!    Ray Leuning with the new radiative transfer scheme   
!    implemented by Y.P. Wang (from Sellers 1986)
!    12/Sept/96 (YPW) correction for mean surface temperature of sunlit
!    and shaded leaves
!    Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)} 
!    
      real Gaussx(5),Gaussw(5)
      real layer1(5),layer2(5)
      real tauL(3),rhoL(3),rhoS(3),Qabs(3,2),Radabv(2),Rnstar(2)
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2),gsleaf(2),QSabs(3,2),Qasoil(2)
      integer ng,nw
      real rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       !Goudriaan

      real rsoil,rlai,raero
      real wsmax,wsmin,WILTPT,FILDCP,wcl(10)
      real gddonset
!    additional arrays to allow output of info for each Layer
      real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
      real GbwcL(5),GswcL(5)
      
! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
!* 5-point
      data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
      data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/

!     soil water conditions
      WILTPT=wsmin/100.
      FILDCP=wsmax/100.
!     reset the vairables
      Rnst1=0.0        !net rad, sunlit
      Rnst2=0.0        !net rad, shaded
      Qcan1=0.0        !vis rad
      Qcan2=0.0
      Rcan1=0.0        !NIR rad
      Rcan2=0.0
      Acan1=0.0        !CO2
      Acan2=0.0
      Ecan1=0.0        !Evap
      Ecan2=0.0
      Hcan1=0.0        !Sens heat
      Hcan2=0.0
      Gbwc1=0.0        !Boundary layer conductance
      Gbwc2=0.0
      Gswc1=0.0        !Canopy conductance
      Gswc2=0.0
      Tleaf1=0.0       !Leaf Temp
      Tleaf2=0.0  
  
!     aerodynamic resistance                                                
      raero=50./wind                           

!    Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam
      
      if(coszen.gt.0) then                                  !check if day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      else
        extKb=100.
      end if

!     Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)
      extkn=extkd                        !N distribution coeff 

!canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
!                                               second; 1=beam, 2=diffuse
      do nw=1,2                                                      !nw:1=VIS, 2=NIR
       scatt(nw)=tauL(nw)+rhoL(nw)                      !scattering coeff
       if((1.-scatt(nw))<0.0)scatt(nw)=0.9999           ! Weng 10/31/2008
       kpr(nw,1)=extKb*sqrt(1.-scatt(nw))               !modified k beam scattered (6.20)
       kpr(nw,2)=extkd*sqrt(1.-scatt(nw))             !modified k diffuse (6.20)
       rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            !canopy reflection black horizontal leaves (6.19)
       rhoc15=2.*xK15*rhoch/(xK15+extkd)                                !canopy reflection (6.21) diffuse
       rhoc45=2.*xK45*rhoch/(xK45+extkd)
       rhoc75=2.*xK75*rhoch/(xK75+extkd)
       rhoc(nw,2)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
       rhoc(nw,1)=2.*extKb/(extKb+extkd)*rhoch                          !canopy reflection (6.21) beam 
       reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))   &                   !effective canopy-soil reflection coeff - beam (6.27)
     &            *exp(-2.*kpr(nw,1)*FLAIT) 
       reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))   &                   !effective canopy-soil reflection coeff - diffuse (6.27)
     &            *exp(-2.*kpr(nw,2)*FLAIT)  
      enddo


!     isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      call Radiso(flait,flait,Qabs,extkd,Tair,eairP,cpair,Patm, &
     &            fbeam,airMa,Rconst,sigma,emleaf,emsoil,       &
     &            emair,Rnstar,grdn)
      TairK=Tair+273.2

!     below      
      do ng=1,5
         flai=gaussx(ng)*FLAIT
!        radiation absorption for visible and near infra-red
         call goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,      &
     &                  scatt,xfang,Qabs) 
!        isothermal net radiation & radiation conductance at canopy top
         call Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,   &
     &               fbeam,airMa,Rconst,sigma,emleaf,emsoil,        &
     &               emair,Rnstar,grdn)
         windUx=wind*exp(-extkU*flai)             !windspeed at depth xi
         scalex=exp(-extkn*flai)                    !scale Vcmx0 & Jmax0
         Vcmxx=Vcmx0*scalex
         eJmxx=eJmx0*scalex
         if(radabv(1).ge.10.0) then                          !check solar Radiation > 10 W/m2
!           leaf stomata-photosynthesis-transpiration model - daytime
            call agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset)       
         else
            call agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
         endif  
         fslt=exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd=1.0-fslt                                !fraction of shaded leaves
         Rnst1=Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
         Rnst2=Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
         RnstL(ng)=Rnst1+Rnst2
!
         Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  !visible
         Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
         QcanL(ng)=Qcan1+Qcan2
!
         Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  !NIR
         Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
         RcanL(ng)=Rcan1+Rcan2
!
         if(Aleaf(1).lt.0.0)Aleaf(1)=0.0      !Weng 2/16/2006
         if(Aleaf(2).lt.0.0)Aleaf(2)=0.0      !Weng 2/16/2006

         Acan1=Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    !amphi/hypostomatous
         Acan2=Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
         AcanL(ng)=Acan1+Acan2

         layer1(ng)=Aleaf(1)
         layer2(ng)=Aleaf(2)

         Ecan1=Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
         Ecan2=Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
         EcanL(ng)=Ecan1+Ecan2
!
         Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
         Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
         HcanL(ng)=Hcan1+Hcan2
!
         Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n
!
         Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
         Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT
      enddo  ! 5 layers

      FLAIT1=(1.0-exp(-extKb*FLAIT))/extkb
      Tleaf1=Tleaf1/FLAIT1
      Tleaf2=Tleaf2/(FLAIT-FLAIT1)

!     Soil surface energy and water fluxes
!    Radiation absorbed by soil
      Rsoilab1=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          !visible
      Rsoilab2=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)        &
     &         +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          !NIR
      Rsoilab1=Rsoilab1*Radabv(1)
      Rsoilab2=Rsoilab2*Radabv(2)
!  
      Tlk1=Tleaf1+273.2
      Tlk2=Tleaf2+273.2
!      temp1=-extkd*FLAIT
      QLair=emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
      QLleaf=emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
     &      +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
      QLleaf=QLleaf*(1.0-exp(-extkd*FLAIT)) 
      QLsoil=emsoil*sigma*(TairK**4)
      Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil

!    Net radiation absorbed by soil
!    the old version of net long-wave radiation absorbed by soils 
!    (with isothermal assumption)
!     Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*FLAIT)         !Longwave
!     Rsoilab3=(1-rhoS(3))*Rsoil3

!    Total radiation absorbed by soil    
      Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3 

!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      fw1=AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.05),1.0)
      Rsoil=30.*exp(0.2/fw1)
      rLAI=exp(FLAIT)
!     latent heat flux into air from soil
!           Eleaf(ileaf)=1.0*
!     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
!     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
     &      (slope+psyc*(rsoil/(raero+rLAI)+1.))
!     sensible heat flux into air from soil
      Hsoil=Rsoilabs-Esoil-G

      return
      end 

!     ****************************************************************************
      subroutine goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,   &
     &                  scatt,xfang,Qabs)
     
!    for spheric leaf angle distribution only
!    compute within canopy radiation (PAR and near infra-red bands)
!    using two-stream approximation (Goudriaan & vanLaar 1994)
!    tauL: leaf transmittance
!    rhoL: leaf reflectance
!    rhoS: soil reflectance
!    sfang XiL function of Ross (1975) - allows for departure from spherical LAD
!         (-1 vertical, +1 horizontal leaves, 0 spherical)
!    FLAI: canopy leaf area index
!    funG: Ross' G function
!    scatB: upscatter parameter for direct beam
!    scatD: upscatter parameter for diffuse
!    albedo: single scattering albedo
!    output:
!    Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
!                       type=1 for sunlit;   =2 for shaded (W/m2)

      real radabv(2)
      real Qabs(3,2),reff(3,2),kpr(3,2),scatt(2)
      xu=coszen                                         !cos zenith angle
      
!     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG=xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam
      
      if(coszen.gt.0) then                                  !check if day or night
        extKb=funG/coszen                                   !beam extinction coeff - black leaves
      else
        extKb=100.
      end if
                       
! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      do nw=1,2
       Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
       Qb0=fbeam*radabv(nw)                                               !beam incident radiation
       Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
     &            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
     &            extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
       Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves 
      end do
      return
      end

!****************************************************************************
      subroutine Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,    &
     &                  fbeam,airMa,Rconst,sigma,emleaf,emsoil,         &
     &                  emair,Rnstar,grdn)
!     output
!     Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
!     23 Dec 1994
!     calculates isothermal net radiation for sunlit and shaded leaves under clear skies
!     implicit real (a-z)
      real Rnstar(2)
      real Qabs(3,2)
      TairK=Tair+273.2

! thermodynamic properties of air
      rhocp=cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K)

! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)       !note eair in Pa
     
! apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)
      tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10

! apparent emissivity from sky plus clouds      
!      emair=emsky+emcloud
! 20/06/96
      emair=emsky

      if(emair.gt.1.0) emair=1.0
      
! net isothermal outgoing longwave radiation per unit leaf area at canopy
! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
      Bn0=sigma*(TairK**4.)
      Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
     &    + exp(-extkd*(flait-flai))*(emsoil-emleaf))
!     isothermal net radiation per unit leaf area for thin layer of sunlit and
!     shaded leaves
      Rnstar(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
      Rnstar(2)=Qabs(1,2)+Qabs(2,2)+Bnxi
!     radiation conductance (m/s) @ flai
      grdn=4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
     &    (exp(-extkd*flai)+exp(-extkd*(flait-flai)))       &
     &    /rhocp
      return
      end
!     ****************************************************************************
      subroutine agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
     &               co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
     &               gsw0,alpha,stom_n,                                 &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,gddonset)

!    implicit real (a-z)
      integer kr1,ileaf
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2), gsleaf(2)
      real Qabs(3,2),Rnstar(2)
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0
!    boundary layer conductance for heat - single sided, forced convection
!    (Monteith 1973, P106 & notes dated 23/12/94)
      if(windUx/wleaf>=0.0)then
          gbHu=0.003*sqrt(windUx/wleaf)    !m/s
      else
          gbHu=0.003 !*sqrt(-windUx/wleaf)
      endif         ! Weng 10/31/2008
!     raero=0.0                        !aerodynamic resistance s/m
      do ileaf=1,2              ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!    ********************************************************************
         kr1=0                     !iteration counter for LE
!        return point for evaporation iteration
         do               !iteration for leaf temperature
!          single-sided boundary layer conductance - free convection (see notes 23/12/94)
           Gras=1.595e8*ABS(Tleaf(ileaf)-Tair)*(wleaf**3.)     !Grashof
           gbHf=0.5*Dheat*(Gras**0.25)/wleaf
           gbH=gbHu+gbHf                         !m/s
           rbH=1./gbH                            !b/l resistance to heat transfer
           rbw=0.93*rbH                          !b/l resistance to water vapour
!          Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
           rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
           rrdn=1./grdn
           Y=1./(1.+ (rbH_L+raero)/rrdn)
!          boundary layer conductance for CO2 - single side only (mol/m2/s)
           gbc=Cmolar*gbH/1.32            !mol/m2/s
           gsc0=gsw0/1.57                 !convert conductance for H2O to that for CO2
           varQc=0.0
           weighR=1.0
           call photosyn(Sps,CO2Ca,CO2Cs,Dleaf,Tlk,Qapar,Gbc,   &   !Qaparx<-Qapar,Gbcx<-Gsc0
     &         theta,a1,Ds0,fwsoil,varQc,weighR,                &
     &         gsc0,alpha,Vcmxx,eJmxx,weighJ,                   &
     &         conKc0,conKo0,Ekc,Eko,o2ci,Rconst,Trefk,         &
     &         Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,       &
     &         Aleafx,Gscx,gddonset)  !outputs
!          choose smaller of Ac, Aq
           Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
!          calculate new values for gsc, cs (Lohammer model)
           co2cs = co2ca-Aleaf(ileaf)/gbc
           co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gscx
!          scale variables
!           gsw=gscx*1.56      !gsw in mol/m2/s, oreginal:gsw=gsc0*1.56,Weng20060215
           gsw=gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
           gswv=gsw/Cmolar                           !gsw in m/s
           rswv=1./gswv
!          calculate evap'n using combination equation with current estimate of gsw
           Eleaf(ileaf)=1.0*                                    &
     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    &   !2* Weng 0215
     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))        

!          calculate sensible heat flux
           Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!          calculate new leaf temperature (K)
           Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
!          calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
           Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
           gbleaf(ileaf)=gbc*1.32*1.075
           gsleaf(ileaf)=gsw
!          compare current and previous leaf temperatures
           if(abs(Tlk1-Tlk).le.0.1) exit ! original is 0.05 C Weng 10/31/2008
!          update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
           Tlk=Tlk1
           Tleaf(ileaf)=Tlk1-273.2
           kr1=kr1+1
           if(kr1 > 500)then
               Tlk=TairK
               exit
           endif
           if(Tlk < 200.)then
                Tlk=TairK
                exit 
           endif                     ! Weng 10/31/2008
!        goto 100                          !solution not found yet
         enddo
! 10  continue
      enddo
      return
      end
!     ****************************************************************************
      subroutine agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,co2ca,    &
     &               wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,            &
     &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,      &
     &               gsw0,alpha,stom_n,                                     &
     &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,                &
     &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,             &
     &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
!    implicit real (a-z)
      integer kr1,ileaf
      real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
      real gbleaf(2), gsleaf(2)
      real Qabs(3,2),Rnstar(2)
!    thermodynamic parameters for air
      TairK=Tair+273.2
      rhocp=cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv=H2oLv0-2.365e3*Tair
      slope=(esat(Tair+0.1)-esat(Tair))/0.1
      psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar=Patm/(Rconst*TairK)
      weighJ=1.0

!     boundary layer conductance for heat - single sided, forced convection
!     (Monteith 1973, P106 & notes dated 23/12/94)
      gbHu=0.003*sqrt(windUx/wleaf)    !m/s
!     raero=0.0                        !aerodynamic resistance s/m

      do ileaf=1,2                  ! loop over sunlit and shaded leaves
!        first estimate of leaf temperature - assume air temp
         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
!        first estimate of deficit at leaf surface - assume Da
         Dleaf=Dair                !Pa
!        first estimate for co2cs
         co2cs=co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1,ileaf)
!        ********************************************************************
         kr1=0                     !iteration counter for LE
         do
!100        continue !    return point for evaporation iteration
!           single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     !Grashof
            gbHf=0.5*Dheat*(Gras**0.25)/wleaf
            gbH=gbHu+gbHf                         !m/s
            rbH=1./gbH                            !b/l resistance to heat transfer
            rbw=0.93*rbH                          !b/l resistance to water vapour
!           Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
            rrdn=1./grdn
            Y=1./(1.+ (rbH_L+raero)/rrdn)
!           boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc=Cmolar*gbH/1.32            !mol/m2/s
            gsc0=gsw0/1.57                        !convert conductance for H2O to that for CO2
            varQc=0.0                  
            weighR=1.0
!           respiration      
            Aleafx=-0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
            gsc=gsc0
!           choose smaller of Ac, Aq
            Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
!           calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca-Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc
!           scale variables
            gsw=gsc*1.56                              !gsw in mol/m2/s
            gswv=gsw/Cmolar                           !gsw in m/s
            rswv=1./gswv
!           calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf)=                                       &
     &      (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/   &
     &      (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
!           calculate sensible heat flux
            Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
!           calculate new leaf temperature (K)
            Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
!           calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf)=gbc*1.32*1.075
            gsleaf(ileaf)=gsw

!          compare current and previous leaf temperatures
            if(abs(Tlk1-Tlk).le.0.1)exit
            if(kr1.gt.500)exit
!           update leaf temperature
            Tlk=Tlk1 
            Tleaf(ileaf)=Tlk1-273.2
            kr1=kr1+1
         enddo                          !solution not found yet
10    continue
      enddo
      return
      end
!     ****************************************************************************
      subroutine ciandA(Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad)
!     calculate coefficients for quadratic equation for ci
      b2 = g0+X*(Gma-Rd)
      b1 = (1.-co2cs*X)*(Gma-Rd)+g0*(Bta-co2cs)-X*(Gma*gammas+Bta*Rd)
      b0 = -(1.-co2cs*X)*(Gma*gammas+Bta*Rd)-g0*Bta*co2cs

      bx=b1*b1-4.*b2*b0
      if(bx.gt.0.0) then 
!       calculate larger root of quadratic
        ciquad = (-b1+sqrt(bx))/(2.*b2)
      endif

      IF(ciquad.lt.0.or.bx.lt.0.) THEN
        Aquad = 0.0
        ciquad = 0.7 * co2Cs
      ELSE
        Aquad = Gma*(ciquad-gammas)/(ciquad+Bta)
      ENDIF
      return
      end

!****************************************************************************
      subroutine goud1(FLAIT,coszen,radabv,fbeam,               &
     &                  Tair,eairP,emair,emsoil,emleaf,sigma,   &
     &                  tauL,rhoL,rhoS,xfang,extkb,extkd,       &
     &                  reffbm,reffdf,extkbm,extkdm,Qcan)
!    use the radiation scheme developed by
!    Goudriaan (1977, Goudriaan and van Larr 1995)
!=================================================================
!    Variable      unit      defintion
!    FLAIT         m2/m2     canopy leaf area index       
!    coszen                  cosine of the zenith angle of the sun
!    radabv(nW)    W/m2      incoming radiation above the canopy
!    fbeam                   beam fraction
!    fdiff                   diffuse fraction
!    funG(=0.5)              Ross's G function
!    extkb                   extinction coefficient for beam PAR
!    extkd                   extinction coefficient for diffuse PAR
!    albedo                  single scattering albedo
!    scatB                   upscattering parameter for beam
!    scatD                   upscattering parameter for diffuse
! ==================================================================
!    all intermediate variables in the calculation correspond
!    to the variables in the Appendix of of Seller (1985) with
!    a prefix of "x".
      integer nW
      real radabv(3)
      real rhocbm(3),rhocdf(3)
      real reffbm(3),reffdf(3),extkbm(3),extkdm(3)
      real tauL(3),rhoL(3),rhoS(3),scatL(3)
      real Qcan(3,2)
!
!     for PAR: using Goudriann approximation to account for scattering
      fdiff=1.0-fbeam
      xu=coszen
      xphi1 = 0.5 -0.633*xfang - 0.33*xfang*xfang
      xphi2 = 0.877 * (1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*xu
      extkb=funG/xu
                       
!     Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180=3.1416/180.
      cozen15=cos(pi180*15)
      cozen45=cos(pi180*45)
      cozen75=cos(pi180*75)
      xK15=xphi1/cozen15+xphi2
      xK45=xphi1/cozen45+xphi2
      xK75=xphi1/cozen75+xphi2
      transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
     &       0.178*exp(-xK75*FLAIT)
      extkd=(-1./FLAIT)*alog(transd)

!     canopy reflection coefficients (Array indices: 1=VIS,  2=NIR
      do nw=1,2                                                         !nw:1=VIS, 2=NIR
         scatL(nw)=tauL(nw)+rhoL(nw)                                    !scattering coeff
         if((1.-scatL(nw))<0.0) scatL(nw)=0.9999                        !Weng 10/31/2008
         extkbm(nw)=extkb*sqrt(1.-scatL(nw))                            !modified k beam scattered (6.20)
         extkdm(nw)=extkd*sqrt(1.-scatL(nw))                            !modified k diffuse (6.20)
         rhoch=(1.-sqrt(1.-scatL(nw)))/(1.+sqrt(1.-scatL(nw)))          !canopy reflection black horizontal leaves (6.19)
         rhoc15=2.*xK15*rhoch/(xK15+extkd)                              !canopy reflection (6.21) diffuse
         rhoc45=2.*xK45*rhoch/(xK45+extkd)
         rhoc75=2.*xK75*rhoch/(xK75+extkd)   
       
         rhocbm(nw)=2.*extkb/(extkb+extkd)*rhoch                        !canopy reflection (6.21) beam 
         rhocdf(nw)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
         reffbm(nw)=rhocbm(nw)+(rhoS(nw)-rhocbm(nw))        &               !effective canopy-soil reflection coeff - beam (6.27)
     &             *exp(-2.*extkbm(nw)*FLAIT)                              
         reffdf(nw)=rhocdf(nw)+(rhoS(nw)-rhocdf(nw))        &            !effective canopy-soil reflection coeff - diffuse (6.27)
     &             *exp(-2.*extkdm(nw)*FLAIT)  

!        by the shaded leaves
         abshdn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &           !absorbed NIR by shaded
     &      *(funE(extkdm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))    ! error found by De Pury
     &      *(funE(extkbm(nw),FLAIT)-funE((extkb+extkbm(nw)),FLAIT))    &
     &      -fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))
!        by the sunlit leaves
         absltn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                       &  !absorbed NIR by sunlit
     &      *funE((extkb+extkdm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-reffbm(nw))*extkbm(nw)                          &
!    &      *funE((extkb+extkdm(nw)),FLAIT)                         ! error found by De Pury
     &      *funE((extkb+extkbm(nw)),FLAIT)                             &
     &      +fbeam*(1.0-scatL(nw))*extkb                                &
     &      *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))

!        scale to real flux 
!        sunlit    
          Qcan(nw,1)=absltn*radabv(nw)
!        shaded
          Qcan(nw,2)=abshdn*radabv(nw)
      enddo
!     
!    calculate the absorbed (iso)thermal radiation
      TairK=Tair+273.2
      
!     apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky=0.642*(eairP/Tairk)**(1./7)      !note eair in Pa

!     apparent emissivity from clouds (Kimball et al 1982)
      ep8z=0.24+2.98e-12*eairP*eairP*exp(3000.0/TairK)
      tau8=amin1(1.0,1-ep8z*(1.4-0.4*ep8z))                !ensure tau8<1
      emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4        !10 from Tcloud = Tair-10 

!     apparent emissivity from sky plus clouds      
!     emair=emsky+emcloud
! 20/06/96
      emair=emsky
      if(emair.gt.1.0) emair=1.0                             

      Bn0=sigma*(TairK**4)
      QLW1=-extkd*emleaf*(1.0-emair)*funE((extkd+extkb),FLAIT)      &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)*exp(-2.0*extkd*FLAIT) &
     &     *funE((extkb-extkd),FLAIT)
      QLW2=-extkd*emleaf*(1.0-emair)*funE(extkd,FLAIT)              &
     &     -extkd*(1.0-emsoil)*(emleaf-emair)                       &
     &     *(exp(-extkd*FLAIT)-exp(-2.0*extkd*FLAIT))/extkd         &
     &     -QLW1
      Qcan(3,1)=QLW1*Bn0
      Qcan(3,2)=QLW2*Bn0
      return
      end

!****************************************************************************
      subroutine photosyn(Sps,CO2Ca,CO2Csx,Dleafx,Tlkx,Qaparx,Gbcx, &
     &         theta,a1,Ds0,fwsoil,varQc,weighR,                    &
     &         g0,alpha,                                            &
     &         Vcmx1,eJmx1,weighJ,conKc0,conKo0,Ekc,Eko,o2ci,       &
     &         Rconst,Trefk,Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,  &
     &         Aleafx,Gscx,gddonset)

!     calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
!     turned on by Weng, 2012-03-13
!     VcmxT = Vjmax(Tlkx,Trefk,Vcmx1,Eavm,Edvm,Rconst,Entrpy)
!     eJmxT = Vjmax(Tlkx,Trefk,eJmx1,Eajm,Edjm,Rconst,Entrpy)
      CO2Csx=AMAX1(CO2Csx,0.6*CO2Ca)
!    check if it is dark - if so calculate respiration and g0 to assign conductance 
      if(Qaparx.le.0.) then                            !night, umol quanta/m2/s
        Aleafx=-0.0089*Vcmx1*exp(0.069*(Tlkx-293.2))   ! original: 0.0089 Weng 3/22/2006
        Gscx=g0
      endif
!     calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
      TminV=gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13
      TmaxV=50.
      ToptV=35.
      
      TminJ=TminV
      TmaxJ=TmaxV
      ToptJ=ToptV 
      
      Tlf=Tlkx-273.2
      VcmxT=VJtemp(Tlf,TminV,TmaxV,ToptV,Vcmx1)
      eJmxT=VJtemp(Tlf,TminJ,TmaxJ,ToptJ,eJmx1)      
!     calculate J, the asymptote for RuBP regeneration rate at given Q
      eJ = weighJ*fJQres(eJmxT,alpha,Qaparx,theta)
!     calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
      conKcT = EnzK(Tlkx,Trefk,conKc0,Rconst,Ekc)
      conKoT = EnzK(Tlkx,Trefk,conKo0,Rconst,Eko)
!     following de Pury 1994, eq 7, make light respiration a fixed proportion of
!     Vcmax
      Rd = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7
      Tdiff=Tlkx-Trefk
      gammas = gam0*(1.+gam1*Tdiff+gam2*Tdiff*Tdiff)       !gamma*
!     gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
      gamma = 0.0
!     ***********************************************************************
!     Analytical solution for ci. This is the ci which satisfies supply and demand
!     functions simultaneously
!     calculate X using Lohammer model, and scale for soil moisture
      a1= 1./(1.-0.7)
      X = a1*fwsoil/((co2csx - gamma)*(1.0 + Dleafx/Ds0))
!     calculate solution for ci when Rubisco activity limits A
      Gma = VcmxT  
      Bta = conKcT*(1.0+ o2ci/conKoT)
      call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci2,Acx)
!     calculate +ve root for ci when RuBP regeneration limits A
      Gma = eJ/4.
      Bta = 2.*gammas
!    calculate coefficients for quadratic equation for ci
      call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci4,Aqx)
!     choose smaller of Ac, Aq
      sps=AMAX1(0.001,sps)                  !Weng, 3/30/2006
      Aleafx = (amin1(Acx,Aqx) - Rd) !*sps     ! Weng 4/4/2006
!      if(Aleafx.lt.0.0) Aleafx=0.0    ! by Weng 3/21/2006
!    calculate new values for gsc, cs (Lohammer model)
      CO2csx = co2ca-Aleafx/Gbcx
      Gscx=g0 + X*Aleafx  ! revised by Weng
      return
      end
!***********************************************************************
      function funeJ(alpha,eJmxT,Qaparx)
      funeJ=alpha*Qaparx*eJmxT/(alpha*Qaparx+2.1*eJmxT)
      return
      end
!****************************************************************************
      real function esat(T)
!     returns saturation vapour pressure in Pa
      esat=610.78*exp(17.27*T/(T+237.3))
      return
      end

!****************************************************************************
      real function evapor(Td,Tw,Patm)
!* returns vapour pressure in Pa from wet & dry bulb temperatures
      gamma = (64.6 + 0.0625*Td)/1.e5
      evapor = esat(Tw)- gamma*(Td-Tw)*Patm
      return
      end

!****************************************************************************
      real function Vjmax(Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
      anum = Vjmax0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      aden = 1. + EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
      Vjmax = anum/aden
      return
      end
!****************************************************************************
      real function funE(extkbd,FLAIT)
      funE=(1.0-exp(-extkbd*FLAIT))/extkbd
      return
      end

!     ****************************************************************************
!     Reed et al (1976, J appl Ecol 13:925) equation for temperature response
!     used for Vcmax and Jmax
      real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
      if(Tlf.lt.TminVJ) Tlf=TminVJ   !constrain leaf temperatures between min and max
      if(Tlf.gt.TmaxVJ) Tlf=TmaxVJ
      pwr=(TmaxVJ-ToptVJ)/(ToptVj-TminVj)
      VJtemp=VJmax0*((Tlf-TminVJ)/(ToptVJ-TminVJ))*     &
     &       ((TmaxVJ-Tlf)/(TmaxVJ-ToptVJ))**pwr 
      return
      end

!     ****************************************************************************
      real function fJQres(eJmx,alpha,Q,theta)
      AX = theta                                 !a term in J fn
      BX = alpha*Q+eJmx                          !b term in J fn
      CX = alpha*Q*eJmx                          !c term in J fn
      if((BX*BX-4.*AX*CX)>=0.0)then
          fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
      else
          fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
      endif

      return
      end

!     *************************************************************************
      real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)

      temp1=(Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk)
!      if (temp1<50.)then
      EnzK = EnzK0*EXP((Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk))
!      else
!      EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
!      endif

      return
      end

!     *************************************************************************
      real function sinbet(doy,lat,pi,timeh)
      real lat
!     sin(bet), bet = elevation angle of sun
!     calculations according to Goudriaan & van Laar 1994 P30
      rad = pi/180.
!     sine and cosine of latitude
      sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
!     sine of maximum declination
      sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
      cosdec=sqrt(1.-sindec*sindec)
!     terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A+B*cos(pi*(timeh-12.)/12.)
      return
      end

!     *************************************************************************
      subroutine yrday(doy,hour,lat,radsol,fbeam)
      real lat
      pi=3.14159256
      pidiv=pi/180.0
      slatx=lat*pidiv
      sindec=-sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
      cosdec=sqrt(1.-sindec*sindec)
      a=sin(slatx)*sindec
      b=cos(slatx)*cosdec
      sinbet=a+b*cos(2*pi*(hour-12.)/24.)
      solext=1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet
      
      tmprat=radsol/solext

      tmpR=0.847-1.61*sinbet+1.04*sinbet*sinbet
      tmpK=(1.47-tmpR)/1.66
      if(tmprat.le.0.22) fdiff=1.0
      if(tmprat.gt.0.22.and.tmprat.le.0.35) then
        fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
      endif
      if(tmprat.gt.0.35.and.tmprat.le.tmpK) then
        fdiff=1.47-1.66*tmprat
      endif
      if(tmprat.ge.tmpK) then
        fdiff=tmpR
      endif
      fbeam=1.0-fdiff
      if(fbeam.lt.0.0) fbeam=0.0
      return
      end


!===============================================================================
!========================================================
! Subroutine 1. Read parameters from file
    subroutine Getparameters(lat,longi,wsmax,wsmin,     &              
    &   LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,           &
    &   SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n,         &
    &   a1,Ds0,Vcmx0,extkU,xfang,alpha,                 &
    &   Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,         &
    &   Tau_Micro,Tau_slowSOM,Tau_Passive,              &
    &   gddonset,Q10,Rl0,Rs0,Rr0,parafile)
    
    implicit none
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmx0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset, Q10,Rl0,Rs0,Rr0
    character(len=50) parafile,commts
    
    parafile=TRIM(parafile)
!   open and read input file for getting climate data
    
    open(10,file=parafile,status='old')
    read(10,11)commts
    read(10,*)lat,longi,wsmax,wsmin
    read(10,11)commts
    read(10,*)LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax    
    read(10,11)commts
    read(10,*)SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    read(10,11)commts
    read(10,*)a1,Ds0,Vcmx0,extkU,xfang,alpha
    read(10,11)commts
    read(10,*)Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,Tau_Micro,Tau_slowSOM,Tau_Passive
    read(10,11)commts
    read(10,*)gddonset,Q10,Rl0,Rs0,Rr0
11  format(a132)
    close(10)
    return
    end

!================================================================  
! Subroutine 1.1 Read estimated parameters      
    subroutine Getparaest(paraestfile,paraest,seq,npara,indexstring)
    implicit none
                
    character(len=50) paraestfile
    integer seq,m,n,istat6
    real paraest(19,40000)
    integer npara
    character(len=250) indexstring

    paraestfile=TRIM(paraestfile)
    open(15,file=paraestfile,status='old',ACTION='read',     &
    &     IOSTAT=istat6)

    read(15,*) npara
    read(15,'(A)') indexstring


    m=0
!   open and read input file for getting climate data
    do
    m=m+1
    read(15,*,IOSTAT=istat6)(paraest(n,m),n=1,(npara+1))
    if(istat6<0)exit
    enddo
    seq=m-1
    close(15)
    return
    end
    
! ====================================================================
! Subroutine 2. Read climatic forcing from file   
    subroutine Getclimate(year_seq,doy_seq,hour_seq,          &
    &   forcing_data,climatefile,lines,yr_length)
    
    implicit none
    integer, parameter :: ilines=150000
    integer, parameter :: iiterms=7
    integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real forcing_data(iiterms,ilines)
    character(len=150) climatefile,commts
    integer m,n,istat1,lines,yr_length

    open(11,file=climatefile,status='old',ACTION='read',     &
    &     IOSTAT=istat1)
!     skip 2 lines of input met data file
      read(11,'(a160)') commts
    m=0  ! to record the lines in a file
    yr_length=0 ! to record years of a dataset
    do    ! read forcing files
        m=m+1
        read(11,*,IOSTAT=istat1)year_seq(m),      &
        &       doy_seq(m),hour_seq(m),           &
        &       (forcing_data(n,m),n=1,iiterms)
        if(istat1<0)exit
    enddo ! end of reading the forcing file
    lines=m-1
    yr_length=(year_seq(lines)-year_seq(1))+1
    close(11)    ! close forcing file
    return
    end
    
!===============================================================================
! Subroutine 3. read observation data from files
    subroutine GetObsData(obs_spruce,std,len1,obsfile1)
         
    implicit none
    real obs_spruce(12,1000),std(12,1000)
    real foliage,fnpp,wood,wnpp,root,rnpp,soilc,phenology
    real foliage_sd,fnpp_sd,wood_sd,wnpp_sd,root_sd,rnpp_sd
    real soilc_sd,pheno_sd,gpp_sd,nee_sd,er_sd
    real gpp,nee,er,sw
    integer days,m,istat2,len1
    integer year,doy,hour,month
    character(len=80) commts,obsfile1

    open(12,file=obsfile1,status='old')
    
!   Observed hourly data
    read(12,901) commts
    m=0
    do
        read (12,*,IOSTAT=istat2)days,gpp,gpp_sd,nee,nee_sd,er,er_sd,   &
        &   foliage,foliage_sd,fnpp,fnpp_sd,wood,wood_sd,wnpp,wnpp_sd,  &
        &   root,root_sd,rnpp,rnpp_sd,soilc,soilc_sd,phenology,pheno_sd
        if(istat2<0)exit
            m=m+1
            obs_spruce(1,m)=real(days)
            obs_spruce(2,m)=gpp
            obs_spruce(3,m)=nee
            obs_spruce(4,m)=er
            obs_spruce(5,m)=foliage
            obs_spruce(6,m)=fnpp
            obs_spruce(7,m)=wood
            obs_spruce(8,m)=wnpp
            obs_spruce(9,m)=root
            obs_spruce(10,m)=rnpp
            obs_spruce(11,m)=soilc
            obs_spruce(12,m)=real(phenology)
            std(1,m)=real(days)
            std(2,m)=gpp_sd
            std(3,m)=nee_sd
            std(4,m)=er_sd
            std(5,m)=foliage_sd
            std(6,m)=fnpp_sd
            std(7,m)=wood_sd
            std(8,m)=wnpp_sd
            std(9,m)=root_sd
            std(10,m)=rnpp_sd
            std(11,m)=soilc_sd
            std(12,m)=real(pheno_sd)
            
            if(days.lt.365*4)then
                len1=m
            endif
    enddo
901 format(a80)  
    close(12)
    return
    end

    
    
!========================================================
! Subroutine: Read Data assimilation check box file
    subroutine GetDAcheckbox(DApar,parmin,parmax,DAparfile)
    
    implicit none
    integer,dimension(35):: DApar
    real,dimension(35):: parmin,parmax
    character(len=50) DAparfile,commts

    DAparfile=TRIM(DAparfile)

    open(15,file=DAparfile,status='old')
    read(15,11)commts
    read(15,*)DApar(1),DApar(2),DApar(3),DApar(4)
    read(15,11)commts
    read(15,*)DApar(5),DApar(6),DApar(7),DApar(8),DApar(9)    
    read(15,11)commts
    read(15,*)DApar(10),DApar(11),DApar(12),DApar(13),DApar(14),DApar(15),DApar(16)
    read(15,11)commts
    read(15,*)DApar(17),DApar(18),DApar(19),DApar(20),DApar(21),DApar(22)
    read(15,11)commts
    read(15,*)DApar(23),DApar(24),DApar(25),DApar(26),DApar(27),DApar(28),DApar(29),DApar(30)
    read(15,11)commts
    read(15,*)DApar(31),DApar(32),DApar(33),DApar(34),DApar(35)
    
    read(15,11)commts
    read(15,*)parmin(1),parmin(2),parmin(3),parmin(4)
    read(15,11)commts
    read(15,*)parmin(5),parmin(6),parmin(7),parmin(8),parmin(9)    
    read(15,11)commts
    read(15,*)parmin(10),parmin(11),parmin(12),parmin(13),parmin(14),parmin(15),parmin(16)
    read(15,11)commts
    read(15,*)parmin(17),parmin(18),parmin(19),parmin(20),parmin(21),parmin(22)
    read(15,11)commts
    read(15,*)parmin(23),parmin(24),parmin(25),parmin(26),parmin(27),parmin(28),parmin(29),parmin(30)
    read(15,11)commts
    read(15,*)parmin(31),parmin(32),parmin(33),parmin(34),parmin(35)
    
    read(15,11)commts
    read(15,*)parmax(1),parmax(2),parmax(3),parmax(4)
    read(15,11)commts
    read(15,*)parmax(5),parmax(6),parmax(7),parmax(8),parmax(9)    
    read(15,11)commts
    read(15,*)parmax(10),parmax(11),parmax(12),parmax(13),parmax(14),parmax(15),parmax(16)
    read(15,11)commts
    read(15,*)parmax(17),parmax(18),parmax(19),parmax(20),parmax(21),parmax(22)
    read(15,11)commts
    read(15,*)parmax(23),parmax(24),parmax(25),parmax(26),parmax(27),parmax(28),parmax(29),parmax(30)
    read(15,11)commts
    read(15,*)parmax(31),parmax(32),parmax(33),parmax(34),parmax(35)
11  format(a132)
    close(15)
    return
    end
    
    
    
! **********************************************************    
    subroutine getCov(gamma,covfile,npara)
    implicit none
    integer npara,i,k
    real gamma(npara,npara)    
    character(len=80) covfile
    
    open(14,file=covfile,status='old')

    do i=1,npara
        read (14,*)(gamma(i,k),k=1,npara)
    enddo  
    return
    end        

      
      
!========================================================================
!    	cost function for observed data
    
    subroutine costFObsNee(Simu_dailyflux,   &
            &   obs_spruce,std,len1,   &
            &   J_last,upgraded)
    
    implicit none
    real Simu_dailyflux(12,10000)
    real obs_spruce(12,1000),std(12,1000)
    integer day
    real J_new,J_last,delta_J,var2
    real J_gpp,J_nee,J_er,J_sw,J_foliage,J_fnpp
    real J_wood,J_wnpp,J_root,J_rnpp,J_soilc,J_pheno
    real tmp1,dObsSim,random_harvest
    integer i,upgraded,len1,len2
    integer j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11
    real r_num
!    compute J_obs
    
    J_gpp=0.0
    J_nee=0.0
    J_er=0.0
    J_foliage=0.0
    J_fnpp=0.0
    J_wood=0.0
    J_wnpp=0.0
    J_root=0.0
    J_rnpp=0.0
    J_soilc=0.0
    J_pheno=0.0
    j1=0
    j2=0
    j3=0
    j4=0
    j5=0
    j6=0
    j7=0
    j8=0
    j9=0
    j10=0
    j11=0
    do i=1,len1
        day=int(obs_spruce(1,i))
        if(obs_spruce(2,i).gt.-999)then
            j1=j1+1
            dObsSim=Simu_dailyflux(1,day)-obs_spruce(2,i)
            J_gpp=J_gpp+(dObsSim*dObsSim)/(2*std(2,i)*std(2,i))
        endif
        if(obs_spruce(3,i).gt.-999)then
            j2=j2+1
            dObsSim=Simu_dailyflux(2,day)-obs_spruce(3,i)
            J_nee=J_nee+(dObsSim*dObsSim)/(2*std(3,i)*std(3,i))
        endif
        if(obs_spruce(4,i).gt.-999)then
            j3=j3+1
            dObsSim=Simu_dailyflux(3,day)-obs_spruce(4,i)
            J_er=J_er+(dObsSim*dObsSim)/(2*std(4,i)*std(4,i))
        endif
        if(obs_spruce(5,i).gt.-999)then
            j4=j4+1
            dObsSim=Simu_dailyflux(4,day)-obs_spruce(5,i)
            J_foliage=J_foliage+(dObsSim*dObsSim)/(2*std(5,i)*std(5,i))
        endif
        if(obs_spruce(6,i).gt.-999)then
            j5=j5+1
            dObsSim=Simu_dailyflux(5,day)-obs_spruce(6,i)
            J_fnpp=J_fnpp+(dObsSim*dObsSim)/(2*std(6,i)*std(6,i))
        endif
        if(obs_spruce(7,i).gt.-999)then
            j6=j6+1
            dObsSim=Simu_dailyflux(6,day)-obs_spruce(7,i)
            J_wood=J_wood+(dObsSim*dObsSim)/(2*std(7,i)*std(7,i))
        endif
        if(obs_spruce(8,i).gt.-999)then
            j7=j7+1
            dObsSim=Simu_dailyflux(7,day)-obs_spruce(8,i)
            J_wnpp=J_wnpp+(dObsSim*dObsSim)/(2*std(8,i)*std(8,i))
        endif
        if(obs_spruce(9,i).gt.-999)then
            j8=j8+1
            dObsSim=Simu_dailyflux(8,day)-obs_spruce(9,i)
            J_root=J_root+(dObsSim*dObsSim)/(2*std(9,i)*std(9,i))
        endif
        if(obs_spruce(10,i).gt.-999)then
            j9=j9+1
            dObsSim=Simu_dailyflux(9,day)-obs_spruce(10,i)
            J_rnpp=J_rnpp+(dObsSim*dObsSim)/(2*std(10,i)*std(10,i))
        endif
        if(obs_spruce(11,i).gt.-999)then
            j10=j10+1
            dObsSim=Simu_dailyflux(10,day)-obs_spruce(11,i)
            J_soilc=J_soilc+(dObsSim*dObsSim)/(2*std(11,i)*std(11,i))
        endif
        if(obs_spruce(12,i).gt.-999)then
            j11=j11+1
            dObsSim=Simu_dailyflux(11,day)-obs_spruce(12,i)
            J_pheno=J_pheno+(dObsSim*dObsSim)/(2*std(12,i)*std(12,i))
        endif
    enddo
 
    J_gpp=J_gpp/real(j1)
    J_nee=J_nee/real(j2)
    J_er=J_er/real(j3)
    J_foliage=J_foliage/real(j4)
    J_fnpp=J_fnpp/real(j5)
    J_wood=J_wood/real(j6)
    J_wnpp=J_wnpp/real(j7)
    J_root=J_root/real(j8)
    J_rnpp=J_rnpp/real(j9)
    J_soilc=J_soilc/real(j10)
    J_pheno=J_pheno/real(j11)
    
    !J_new=J_gpp+J_er+J_foliage+J_wood   
    !J_new = J_gpp+J_er+J_foliage+J_wood+J_fnpp+J_wnpp
    !J_new = J_gpp+J_er+J_foliage+J_wood+J_root+J_soilc
    J_new=J_gpp+J_er+J_foliage+J_wood+J_fnpp+J_wnpp+J_root+J_rnpp+J_soilc+J_pheno 
    !J_new=J_soilc 
    delta_J=J_new-J_last           !    delta_J=(J_new-J_last)/J_last
    
    CALL random_number(r_num)

!     delta_J=-1     !accept all samples, No data
    if(ISNAN(J_new))then
        write(*,*)'NaN return, upgraded', upgraded 
        return
    endif
    if(AMIN1(1.0,exp(-delta_J)).gt.r_num)then
        upgraded=upgraded+1
        J_last=J_new    
    endif
    write(*,*) upgraded,delta_J,J_last
    return
    end


!******************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Square root of a matrix							  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine racine_mat(M, Mrac,npara)

    integer npara,i
    real M(npara,npara),Mrac(npara,npara)
    real valpr(npara),vectpr(npara,npara)
    Mrac=0.
    call jacobi(M,npara,npara,valpr,vectpr,nrot)
    do i=1,npara
	if(valpr(i).ge.0.) then
            Mrac(i,i)=sqrt(valpr(i))
	else
            print*, 'WARNING!!! Square root of the matrix is undefined.'
            print*, ' A negative eigenvalue has been set to zero - results may be wrong'
            Mrac=M
            return
	endif
    enddo
    Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

end subroutine racine_mat      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extraction of the eigenvalues and the eigenvectors !!
!! of a matrix (Numerical Recipes)					  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER :: n,np,nrot
REAL :: a(np,np),d(np),v(np,np)
INTEGER, PARAMETER :: NMAX=500
INTEGER :: i,ip,iq,j
REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      
do ip=1,n
	do iq=1,n
		v(ip,iq)=0.
	end do
	v(ip,ip)=1.
end do

do ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.
end do

nrot=0
do i=1,50
	sm=0.
	do ip=1,n-1
		do iq=ip+1,n
			sm=sm+abs(a(ip,iq))
		end do
	end do
	if(sm.eq.0.)return
	if(i.lt.4)then
		tresh=0.2*sm/n**2
	else
		tresh=0.
	endif
	do ip=1,n-1
		do iq=ip+1,n
			g=100.*abs(a(ip,iq))
			if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
				a(ip,iq)=0.
			else if(abs(a(ip,iq)).gt.tresh)then
				h=d(iq)-d(ip)
				if(abs(h)+g.eq.abs(h))then
					t=a(ip,iq)/h
				else
					theta=0.5*h/a(ip,iq)
					t=1./(abs(theta)+sqrt(1.+theta**2))
					if(theta.lt.0.) then
						t=-t
					endif
				endif
				c=1./sqrt(1+t**2)
				s=t*c
				tau=s/(1.+c)
				h=t*a(ip,iq)
				z(ip)=z(ip)-h
				z(iq)=z(iq)+h
				d(ip)=d(ip)-h
				d(iq)=d(iq)+h
				a(ip,iq)=0.
				do j=1,ip-1
					g=a(j,ip)
					h=a(j,iq)
					a(j,ip)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=ip+1,iq-1
					g=a(ip,j)
					h=a(j,iq)
					a(ip,j)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=iq+1,n
					g=a(ip,j)
					h=a(iq,j)
					a(ip,j)=g-s*(h+g*tau)
					a(iq,j)=h+s*(g-h*tau)
				end do
				do j=1,n
					g=v(j,ip)
					h=v(j,iq)
					v(j,ip)=g-s*(h+g*tau)
					v(j,iq)=h+s*(g-h*tau)
				end do
				nrot=nrot+1
			endif
		end do
	end do
	do ip=1,n
		b(ip)=b(ip)+z(ip)
		d(ip)=b(ip)
		z(ip)=0.
	end do
end do
print*, 'too many iterations in jacobi' 
return
END subroutine jacobi


!===================================================
!       generate new coefficents
        subroutine coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        
        integer npara
        real coefac(npara),coefmax(npara),coefmin(npara),coef(npara)
        real r,coefmid,random_harvest
        integer i
        real search_length
        do i=1,npara
999         continue
            CALL random_number(random_harvest)
            r=random_harvest-0.5
            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length
            if(coef(i).gt.coefmax(i).or.coef(i).lt.coefmin(i))goto 999
        enddo
        return
        end
!============
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random vector from a multivariate  !!
!! normal distribution with mean zero and covariance  !!
!! matrix gamma.									  !!
!! Beware!!! In order to improve the speed of the	  !!
!! algorithms, the subroutine use the Square root	  !!
!! matrix of gamma									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gengaussvect(gamma_racine,xold,xnew,npara)

integer npara
real gamma_racine(npara,npara)
real x(npara),xold(npara),xnew(npara)

do i=1,npara
    x(i)=rangauss(25)
enddo

x = matmul(gamma_racine, x)
xnew = xold + x
end subroutine gengaussvect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random number from a standard	  !!
!! normal distribution. (Numerical Recipes)           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rangauss(idum)


integer idum
real v1, v2, r, fac, gset
real r_num

data iset/0/
if(iset==0) then
1	CALL random_number(r_num)
        v1=2.*r_num-1
        CALL random_number(r_num)
	v2=2.*r_num-1
	r=(v1)**2+(v2)**2
	if(r>=1) go to 1
	fac=sqrt(-2.*log(r)/r)
	gset=v1*fac
	rangauss=v2*fac
	iset=1
else
	rangauss=gset
	iset=0
end if

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variance matrix of a matrix of data				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine varcov(tab,varcovar,npara,ncov)

integer npara,ncov
real tab(ncov,npara),tab2(ncov,npara)
real varcovar(npara,npara)

call centre(tab,tab2,npara,ncov)

varcovar = matmul(transpose(tab2), tab2)*(1./real(ncov))

end subroutine varcov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the centered matrix, ie. the matrix minus  !!
!! the column means									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine centre(mat,mat_out,npara,ncov)

    integer npara,ncov
    real mat(ncov,npara),mat_out(ncov,npara)
    real mean

do i=1,npara
    mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
enddo

end subroutine centre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mean of a vector									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function mean(tab,ncov)
    integer ncov
real tab(ncov)
real mean,mean_tt
mean_tt=0.
do i=1,ncov	
mean_tt=mean_tt+tab(i)/real(ncov)
enddo
mean=mean_tt
End Function

