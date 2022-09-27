       program TECOS
       implicit none
       integer, parameter :: No_site=1     
       integer, parameter :: iiterms=15     
       integer, parameter :: iterms_NEE=13 
       integer, parameter :: ilines=1402560  
       integer, parameter :: dlines=58440  
       real, parameter:: times_storage_use=600  
       integer dayCount 
       integer nonValidCount 
       integer nValid  
       real sValid,varValue,paramspec,slope,sTotal  
       real sumObs,sumSim,sumObsObs,sumSimObs  
       real NEE_diff,bestSum,NEE_dif_sum 
       integer lines,idays
       real inputstep,step_NEE
       integer,dimension(ilines):: year_data,year_NEE
       real,dimension(ilines) :: doy_data,hour_data
       real,dimension(ilines) :: doy_NEE,hour_NEE
       character(len=10) s_dlines
       character(len=10) s_ilines
       real NEE_obs(ilines) 
       real neeSim_daily,neeObs_daily
       real neeDif_day,neeDif_year,neeDif_all 
       real,dimension(dlines)::neeSim_days,neeObs_days 
       real input_data(iiterms,ilines)
       real reader_NEE(iterms_NEE)

       real lat,longi,rdepth,LAIMAX,LAIMIN
       real wsmax,wsmin,co2ca
       real tau_R1,tau_R2,tau_R3 
       real tau_S1,tau_S2,tau_S3 
       real tau_L,tau_W,tau_R
       real tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass

       real Q_leaf,Q_wood,Q_root
       real Q_fine,Q_coarse,Q_Micr,Q_Slow,Q_Pass
       real Rh_f,Rh_c,Rh_Micr,Rh_Slow,Rh_Pass

       real TminV,TmaxV,ToptV,Tcold,Gamma_Wmax,Gamma_Tmax

       character(len=6) vegtype
       character(len=6) site

       real WILTPT,FILDCP
       real fwsoil,topfws,wscontent,omega,omega_s
       real WaterR,WaterS(10),SatFracL(10)

       real NSC,NSCmin,NSCmax,add       
       real Growth,Groot,Gshoot,GRmax   
       real St,Sw,Ss,Sn,Srs,Sps,fnsc,Weight 
       real Twsoil(7),Tavg,Tcur

       real evap,transp,ET,G,Qh,Qle 
       real wind,windu0,eairp,esat,wethr,rnet
       real LWdown,Pa_air,Qair  
       real gpp,NPP,NEE,gpp_t,evap_t,transp_t
       real,dimension(3):: tauL,rhoL,rhoS,reffbm
       real,dimension(3):: reffdf,extkbm,extkdm
       real,dimension(2):: Radabv,Acan,Ecan,Hcan
       real,dimension(2):: Tcan,Gbwcan,Gswcan
       real Qcan(3,2),Qcan0(3)

       real stom_n,a1,Ds0,Vcmx0,extkU,xfang,alpha
       real pi,emleaf,emsoil
       real Rconst,sigma,cpair,Patm,Trefk,H2OLv0
       real airMa,H2OMw,chi,Dheat
       real wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci
       real Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2

       real,dimension(5):: RnStL,QcanL,RcanL,AcanL,EcanL,HcanL
       real,dimension(5):: GbwcL,GswcL,hG,hIL
       real,dimension(5):: Gaussx,Gaussw,Gaussw_cum 

       real,dimension(32):: randnums
       real,dimension(1,45):: params
       real S_w_min,Q10_h 

       real LAI,bmroot,bmstem,bmleaf,bmplant
       real SLA,L_fall,L_add,litter,seeds
       real GDDonset,GDD5,accumulation,storage,stor_use,store
       real RaL,RaS,RaR  
       real alpha_L,alpha_W,alpha_R 
       real RaLeaf,RaStem,RaRoot,Rsoil,Rauto,Rhetero,Rtotal 
       real gpp_yr,NPP_yr,NEE_yr,RaL_yr,RaR_yr,RaS_yr,Rh_yr
       real NPPL_yr,NPPR_yr,NPPS_yr,NPP_L,NPP_R,NPP_W
       real Rootmax,Stemmax,SapS,SapR,StemSap,RootSap
       REAL ws,wdepth

       real Ta,Tair,Tavg72,Ts
       real doy,hour,Tsoil,Dair,Rh,rain,radsol,rain_t

       real CO2air_d_avg,LWdown_d_avg,SWdown_d_avg,Psurf_d_avg
       real Qair_d_avg,Rain_d_avg,Tair_d_avg,Wind_d_avg


       real Q_root1,Q_root2,Q_root3 
       real Q_soil,Q_soil1,Q_soil2,Q_soil3 
       real PlantC

       real TEVAP,AEVAP,evap_yr,transp_yr
       real,dimension(10):: thksl,wupl,evapl,wcl,FRLEN
       real runoff,Trunoff,runoff_yr,rain_yr
       real wsc(10),ws1,ws2,dws,net_dws
       real gc_wengL,gp_wengL,gc_wengD,gp_wengD
       real gcgp_D,gcgp_yr,gcgp
       real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop
       real Hcanop,Hcanop_d
       real Raplant,Glmax,Gsmax,Rh_d,ET_d

       real Sc_co2,Sc_T,Sc_prcp,CO2

       real tranformer
       real NACP_L,NACP_W,NACP_R,NACP_F
       real NACP_C,NACP_S1,NACP_S2,NACP_S3
       real NACP_GPP,NACP_NPP,NACP_NEE
       real NACP_Ra,NACP_Rh,NACP_Rt   
       real NACP_ET,NACP_TR,NACP_Run
       real NA_L,NA_W,NA_R,NA_can,NA_BM
       real NA_C,NA_F,NA_Micr,NA_Slow,NA_Pass,NA_S
       real NA_GPP,NA_NPP,NA_NEE
       real NA_Ra,NA_Rh,NA_Rt
       real NA_ET,NA_TR,NA_Run
       real NA_WC_soil,NA_WC_root,NA_wcLayer(10) 

       real NEE_annual,Cumol2gram
       real NEE_annual_array(30)
       integer year_array(30),year_obs

       integer jrain,num_gcgp,W_flag(7)
       integer onset,duration,offset,dormancy  
       integer year,yr,days,i,j,k,m,n,irunmean,writer
       integer lines_NEE,yr_NEE
       integer istat1,istat2,istat3,istat4
       integer dtimes,yr_data
       integer num_scen,isite
       integer idoy,ihour,ileaf,num

       character(len=150) climfile,NEE_file,out_yr,out_d
       character(len=150) parafile,logfile       
       character(len=400) Variables,Units       
       character(len=80) commts


       call getarg(1,parafile)
       open(10,file=parafile,status='old',ACTION='read',IOSTAT=istat1)
       if(istat1==0)then
       else
         close(10)
         goto 9999
       endif

       read(10,11)commts
       do isite=1,No_site  
         read(10,*,IOSTAT=istat1)site,vegtype,climfile,NEE_file,out_d,
     &     lat,longi,wsmax,wsmin,gddonset,LAIMAX,LAIMIN,rdepth,
     &     Rootmax,Stemmax,SapR,SapS,SLA,GLmax,GRmax,Gsmax,
     &     stom_n,a1,Ds0,Vcmx0,extkU,xfang,alpha,co2ca,
     &     tau_L,tau_W,tau_R,
     &     tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass             
         if(istat1 /= 0)then
           goto 999
         endif
       close(10) 

       call getarg(2,out_d)



         tau_L =tau_L *8760.                          
         tau_W =tau_W *8760.
         tau_R =tau_R *8760.
         tau_F =tau_F *8760.
         tau_C =tau_C *8760.
         tau_Micr=tau_Micr*8760.
         tau_Slow=tau_Slow*8760.
         tau_Pass=tau_Pass*8760.

         GLmax=GLmax/24.
         GRmax=GRmax/24.
         Gsmax=GSmax/24.

      dayCount = 1
         nonValidCount=0
         sumObs=0.0	
         sumSim=0.0
         sumObsObs=0.0
         sumSimObs=0.0

         co2=co2ca
         site=trim(site)
         vegtype=trim(vegtype)
         climfile=trim(climfile)
         open(11,file=climfile,status='old',
     &     ACTION='read',IOSTAT=istat2)
         if(istat2==0)then
         else
           close(11)
           goto 999
         endif

         read(11,'(a160)') commts
         read(11,'(a160)') commts
11       format(a160)
         m=0  
         yr_data=0 

         do    
           m=m+1
           read(11,*,IOSTAT=istat3)year_data(m),
     &       doy_data(m),hour_data(m),
     &       (input_data(n,m),n=1,iiterms)
           hour_data(m)=hour_data(m)+1
           if(istat3<0)exit
         enddo 
         close(11)    

         lines=m-1
         yr_data=(year_data(lines)-year_data(1))+1
         inputstep=hour_data(2)-hour_data(1)
         if (inputstep==1.0)then
         else
           goto 999
         endif



         open(12,file=NEE_file,status='old',ACTION='read',
     &     IOSTAT=istat2)
         if(istat2==0)then
           read(12,11) commts
           m=0
           do 
             m=m+1
             read(12,*,IOSTAT=istat3)year_NEE(m),doy_NEE(m),
     &         hour_NEE(m),(reader_NEE(n),n=1,iterms_NEE)
             if(istat3<0)exit
             if(reader_NEE(11)==-999.0)reader_NEE(11)=0
             NEE_obs(m)=reader_NEE(11)*12./1000./1000000.   
           enddo
         endif
         lines_NEE=m-1
         yr_NEE=(year_NEE(lines_NEE)-year_NEE(1))+1
         step_NEE=hour_NEE(2)-hour_NEE(1)
         close(12)

         Cumol2gram=3600.*Step_NEE*1000
         NEE_annual=0.0
         m=0
         year_obs=year_NEE(1)
         do i=1,lines_NEE
           if(year_obs.eq.year_NEE(i)) then
             if(NEE_obs(i)<-500.*12./1000./1000000)then
               NEE_obs(i)=0.0
               nonValidCount=nonValidCount+1
             else
               NEE_annual=NEE_annual+NEE_obs(i)*Cumol2gram
             end if
           else
             m=m+1
             NEE_annual_array(m)=NEE_annual
             year_array(m)=year_obs
             NEE_annual=0
             year_obs=year_NEE(i)
             NEE_annual=NEE_annual+NEE_obs(i)*Cumol2gram
           endif
         enddo

         open(21,file=out_d)
         write(21,'(56(a15,1X))')'year','doy','hour',
     &     'carbon_Leaf','carbon_Wood','carbon_Root','carbon_Flitter',
     &     'carbon_Clitter','SOM_Miro','SOM_SLOW','SOM_Pass',
     &     'carbon_canopy','carbon_biomass','carbon_soil',
     &     'gpp','npp','nee_simulate','nee_observe',
     &     'resp_auto','resp_hetero','resp_tot',
     &     'ET','Transpiration','Runoff','LatentHeat','LAI','RootMoist',
     &     'SoilWater','soilwater_1','soilwater_2','soilwater_3',
     &     'soilwater_4','soilwater_5','soilwater_6','soilwater_7',
     &     'soilwater_8','soilwater_9','soilwater_10',
     &     'satfrac_1','satfrac_2','satfrac_3','satfrac_4','satfrac_5',
     &     'satfrac_6','satfrac_7','satfrac_8','satfrac_9','satfrac_10',
     &     'CO2concentration','AirSpecificHumudity','rain_div_3600',
     &     'scale_sw','TairPlus273_15'




         thksl(1)=10.0   
         thksl(2)=10.0   
         thksl(3)=10.0   
         thksl(4)=10.0   
         thksl(5)=10.0   
         thksl(6)=20.0   
         thksl(7)=20.0   
         thksl(8)=20.0   
         thksl(9)=20.0   
         thksl(10)=20.0  

         FRLEN(1)=0.30
         FRLEN(2)=0.20
         FRLEN(3)=0.15
         FRLEN(4)=0.15
         FRLEN(5)=0.1
         FRLEN(6)=0.05
         FRLEN(7)=0.05
         FRLEN(8)=0.0
         FRLEN(9)=0.0
         FRLEN(10)=0.0

      call consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,
     &   Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,
     &   wleaf,gsw0,Vcmx0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)

         OPEN(110, FILE="initial_opt.txt")
         READ(110,*)PARAMS
         close(110)






         wsmax=PARAMS(1,1)
         wsmin=PARAMS(1,2)
         gddonset=PARAMS(1,3)
         LAIMAX=PARAMS(1,4)
         LAIMIN=PARAMS(1,5)
         rdepth=PARAMS(1,6)
         Rootmax=PARAMS(1,7)
         Stemmax=PARAMS(1,8)
         SapR=PARAMS(1,9)
         SapS=PARAMS(1,10)
         SLA=PARAMS(1,11)
         GLmax=PARAMS(1,12)/24.
         GRmax=PARAMS(1,13)/24.
         Gsmax=PARAMS(1,14)/24.
         a1=PARAMS(1,15)
         Ds0=PARAMS(1,16)
         Vcmx0=PARAMS(1,17)
         alpha=PARAMS(1,18)
         tau_L=PARAMS(1,19)*8760.
         tau_W=PARAMS(1,20)*8760.
         tau_R=PARAMS(1,21)*8760.
         tau_F=PARAMS(1,22)*8760.
         tau_C=PARAMS(1,23)*8760.
         tau_Micr=PARAMS(1,24)*8760.
         tau_Slow=PARAMS(1,25)*8760.
         tau_Pass=PARAMS(1,26)*8760.
         TminV=PARAMS(1,27)
         TmaxV=PARAMS(1,28)
         ToptV=PARAMS(1,29)
         Tcold=PARAMS(1,30)
         Gamma_Wmax=PARAMS(1,31)
         Gamma_Tmax=PARAMS(1,32)


         nsc=PARAMS(1,33)
         Q_leaf=PARAMS(1,34)
         Q_wood=PARAMS(1,35)
         Q_root1=PARAMS(1,36)
         Q_root2=PARAMS(1,37)
         Q_root3=PARAMS(1,38)
         Q_coarse=PARAMS(1,39)
         Q_fine=PARAMS(1,40)
         Q_micr=PARAMS(1,41)
         Q_slow=PARAMS(1,42)
         Q_pass=PARAMS(1,43)
         S_w_min=PARAMS(1,44)
         Q10_h=PARAMS(1,45)



         WILTPT=wsmin/100.0
         FILDCP=wsmax/100.0

         wscontent=WILTPT
         fwsoil=1.0
         topfws=1.0
         omega=1.0
         do i=1,10
           wcl(i)=FILDCP
         enddo
         Storage=32.09           
         stor_use=Storage/times_storage_use
         onset=0
         duration=0
         offset=0
         dormancy=1 

         LAI=LAIMIN
         bmstem=Q_wood/0.45
         bmroot=(Q_root1+Q_root2+Q_root3)/0.45
         bmleaf=Q_leaf/0.45
         bmplant=bmstem+bmroot+bmleaf



         writer=yr_data*1  
         num=0
         m=1
         n=1
         Tavg72=5.0
         neeDif_all=0.0 
         do yr=1,writer+yr_data  






           if(mod(year_data(m), 4)==0) then
               if(mod(year_data(m),100)==0) then
                   if(mod(year_data(m),400)==0) then
                       idays=366
                   else
                       idays=365
                   endif
               else
                   idays=366
               endif
           else
               idays=365
           endif

           GDD5=0.0
           onset=0
           gpp_yr=0.0
           NPP_yr=0.0
           Rh_yr =0.0
           NEE_yr=0.0
           neeDif_year=0.0
           do days=1,idays 
             StemSap=AMIN1(Stemmax,SapS*bmStem)
             RootSap=AMIN1(Rootmax,SapR*bmRoot)
             NSCmin=5. 
             NSCmax=0.05*(StemSap+RootSap+Q_leaf)
             if(Ta.gt.0.0)GDD5=GDD5+Ta-0.0 

             gpp_t   =0.0   
             transp_t=0.0   
             Hcanop_d=0.0   
             evap_t  =0.0   
             ta=0.0         
             Ts=0.0         
             rain_t=0.0     
             Trunoff=0.0    
             RaL=0.0
             RaS=0.0
             RaR=0.0
             Rauto=0.0
             neeDif_day=0.0 
             neeSim_daily=0.0 
             neeObs_daily=0.0 
             dtimes=24 
             do i=1,dtimes

               if(m > lines)then 
                 m=1
                 n=1
               endif
               year=year_data(m)
               doy=doy_data(m)
               hour =hour_data(m)
               Tair=input_data(1,m)-273.15   
               rain=input_data(7,m)*3600.    
               radsol=input_data(11,m)       
               Qair=input_data(3,m)          
               wind=ABS(input_data(5,m))     
               co2ca=input_data(15,m) *1.0E-6 
               Pa_air=101325.0   
               LWdown=input_data(13,m)
               if(inputstep.eq.1.0)then
                 m=m+1
               else
                 rain=(input_data(7,m)+input_data(7,m+1))*3600.
                 m=m+2
               endif
               n=n+1
201            format(I4,(1x,f4.0),6(1x,f8.4))
               Tsoil=Tair*0.8    
               windU0=Min(1.0,wind) 
               if(radsol.eq.0.0) radsol=0.01


               eairP=Qair/(Qair+0.62198)*Pa_air         
               Dair=Max(50.01,(esat(Tair)-eairP))       
               RH=eairP/(eairP+Dair)*100.0
               wethr=1
               Rnet=0.8*radsol
               if(radsol.gt.10.0) then
                 G=-25.0
               else
                 G=20.5
               endif
               Esoil=0.05*radsol
               if(radsol.LE.10.0) Esoil=0.5*G
               Hcrop=0.1  
               Ecstot=0.1 
               Anet=0.1 
               DepH2O=0.2

               ta= ta + tair/24.0             
               Ts=Ts+Tsoil/24.0

               if(NSC.le.NSCmin)fnsc=0.0
               if(NSC.ge.NSCmax)fnsc=1.0
               if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                 fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
               endif
               call canopy(gpp,evap,transp,Acanop,Hcanop,   
     &           fwsoil,topfws,wscontent,                    
     &           LAI,Sps,
     &           doy,hour,radsol,tair,dair,eairP,            
     &           windU0,rain,wethr,
     &           Rnet,G,Esoil,Hcrop,Ecstot,Anet,
     &           Tsoil,DepH2O,
     &           wsmax,wsmin,                                
     &           lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,
     &           stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,
     &           Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,
     &           H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,
     &           conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,
     &           Edjm,Entrpy,gam0,gam1,gam2,TminV,TmaxV,ToptV)

               gpp=Amax1(0.0,gpp)
               call respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,
     &           LAI,SLA,bmstem,bmroot,bmleaf,
     &           StemSap,RootSap,NSC,fnsc,
     &           RaLeaf,RaStem,RaRoot,Rauto)

               call soilwater(wsmax,wsmin,rdepth,FRLEN,
     &           rain,tair,transp,wcl,tsoil,Rh,thksl,LAI,   
     &           evap,runoff,wscontent,fwsoil,topfws,  
     &           omega,omega_S,WaterR,WaterS,SatFracL)  
               omega=omega_S
               ET=evap+transp
               NA_WC_soil=wscontent
               NA_WC_root=WaterR
               NA_wcLayer=SatFracL

               call plantgrowth(Tavg72,Tavg72,omega,GLmax,GRmax,
     &           GSmax,LAI,LAIMAX,LAIMIN,SLA,Tau_L,
     &           bmleaf,bmroot,bmstem,bmplant,
     &           Rootmax,Stemmax,SapS,
     &           StemSap,RootSap,Storage,GDD5,
     &           stor_use,onset,accumulation,gddonset,
     &           Sps,NSC,fnsc,NSCmin,NSCmax,
     &           store,add,L_fall,Tcold,Gamma_Wmax,Gamma_Tmax,
     &           NPP,alpha_L,alpha_W,alpha_R)

               NSC=NSC+GPP-Rauto-(NPP-add)-store
               write(*,*)GPP,Rauto,npp,add,store


               call TCS(Tair,Tsoil,omega,
     &           NPP,alpha_L,alpha_W,alpha_R,
     &           L_fall,tau_L,tau_W,tau_R,
     &           tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass,
     &           Q_leaf,Q_wood,Q_root,
     &           Q_fine,Q_coarse,Q_Micr,Q_Slow,Q_Pass,
     &           Rh_f,Rh_c,Rh_Micr,Rh_Slow,Rh_Pass,S_w_min,Q10_h)

               Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
               NEE=Rauto+Rhetero - GPP
               Q_soil=Q_Micr + Q_Slow + Q_Pass
               Q_root=Q_root
               bmroot=Q_root/0.45
               bmleaf=Q_leaf/0.45
               bmstem=Q_wood/0.45
               bmplant=bmleaf+bmroot+bmstem
               LAI=bmleaf*SLA

               if((yr.gt.writer).and.(yr.le.(writer+yr_data)))then
                 NA_L=Q_leaf/1000. 
                 NA_W=Q_wood/1000. 
                 NA_R=Q_root/1000. 
                 NA_F=Q_fine/1000. 
                 NA_C=Q_coarse/1000. 
                 NA_Micr=Q_Micr/1000. 
                 NA_slow=Q_Slow/1000. 
                 NA_pass=Q_Pass/1000. 
                 NA_can=NA_L+NA_W*0.7
                 NA_BM=NA_L+NA_W+NA_R
                 NA_S=NA_Micr+NA_slow+NA_pass 
                 NA_GPP=GPP/(1000.*3600.) 
                 NA_NPP=NPP/(1000.*3600.) 
                 NA_NEE=NEE/(1000.*3600.) 
                 NA_Ra=Rauto/(1000.*3600.) 
                 NA_Rh=Rhetero/(1000.*3600.) 
                 NA_Rt=NA_Ra+NA_Rh
                 NA_ET=ET            
                 NA_TR=transp
                 NA_Run=runoff
                 Qle=ET*((2.501-0.00236*Tair)*1000000.0)/3600. 
                 if(NEE_obs(n)>-500.*12./1000./1000000) then 
                   NEE_diff = (NA_NEE-NEE_obs(n))*100000000/12. 
                 else
                   NEE_diff = 0.0                 
                 endif
                   write(21,121)year,doy,hour,
     &               NA_L,NA_W,NA_R,NA_F,NA_C,
     &               NA_Micr,NA_slow,NA_pass,
     &               NA_can,NA_BM,NA_S,                        
     &               NA_GPP,NA_NPP,NA_NEE,NEE_obs(n),
     &               NA_Ra,NA_Rh,NA_Rt,   
     &               NA_ET,NA_TR,NA_Run,Qle,
     &               LAI,WaterR,wscontent,
     &               (WaterS(k),k=1,10),(SatFracL(k),k=1,10),
     &               CO2ca,Qair,rain/3600.,radsol,
     &               Tair+273.15            
               endif

               gpp_t=gpp_t + gpp*(24./dtimes)
               transp_t=transp_t + transp*(24./dtimes)
               evap_t=evap_t + evap*(24./dtimes)
               Hcanop_d=Hcanop_d+Hcanop/(24./dtimes)
               Trunoff=Trunoff+runoff

               gpp_yr=gpp_yr+gpp
               NPP_yr=NPP_yr+NPP
               Rh_yr =Rh_yr +Rhetero
               NEE_yr=NEE_yr+NEE
               if((yr.gt.writer).and.(yr.le.(writer+yr_data)))then
                 neeDif_day= neeDif_day + NEE_diff*NEE_diff 
                 neeSim_daily=neeSim_daily+NA_NEE*1000000000./12. 
                 neeObs_daily=neeObs_daily+NEE_obs(n)*1000000000./12. 
               endif
                                   
             enddo              
             if((yr.gt.writer).and.(yr.le.(writer+yr_data)))then
               neeObs_days(dayCount) = neeObs_daily
               neeSim_days(dayCount) = neeSim_daily
               sumObs=sumObs+neeObs_daily
               sumSim=sumSim+neeSim_daily
               sumObsObs=sumObsObs+neeObs_daily*neeObs_daily
               sumSimObs=sumSimObs+neeSim_daily*neeObs_daily
           dayCount = dayCount + 1
               neeDif_year=neeDif_year + neeDif_day 
         endif

             Tavg72=ta
121          format(I6,1X,f15.0,1X,f15.2,1X,120(E15.8,1X))
           enddo                         
           storage=accumulation
           stor_use=Storage/times_storage_use
           accumulation=0.0
           onset=0
           neeDif_all=neeDif_all + neeDif_year 
         enddo            
         close(21)
999      continue
       enddo   

9999   continue

      paramspec=0.0
      sTotal=0.0
      nValid = lines_NEE-nonValidCount
      varValue = neeDif_all/nValid
      sValid = sqrt(varValue)



       end



      subroutine consts(pi,tauL,rhoL,rhoS,emleaf,emsoil,
     &   Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat,
     &   wleaf,gsw0,Vcmx0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2)
     
      real tauL(3), rhoL(3), rhoS(3)
      pi = 3.1415926

      tauL(1)=0.1                  
      rhoL(1)=0.1                  
      rhoS(1)=0.1                  
      tauL(2)=0.425                
      rhoL(2)=0.425                
      rhoS(2)=0.3                  
      tauL(3)=0.00                 
      rhoL(3)=0.00                 
      rhoS(3)=0.00                 
      emleaf=0.96
      emsoil=0.94
      Rconst=8.314                 
      sigma=5.67e-8                
      cpair=1010.                  
      Patm=1.e5                    
      Trefk=293.2                  
      H2OLv0=2.501e6               
      AirMa=29.e-3                 
      H2OMw=18.e-3                 
      chi=0.93                     
      Dheat=21.5e-6                

      gsw0 = 1.0e-2                
      eJmx0 = Vcmx0*2.7            
      theta = 0.9
      wleaf=0.01                   


      conKc0 = 302.e-6                
      conKo0 = 256.e-3                
      Ekc = 59430.                    
      Eko = 36000.                    
      o2ci= 210.e-3                   


      Eavm = 116300.               
      Edvm = 202900.               
      Eajm = 79500.                
      Edjm = 201000.               
      Entrpy = 650.                


      gam0 = 28.0e-6               
      gam1 = .0509
      gam2 = .0010
      return
      end




       subroutine canopy(gpp,evap,transp,Acanop,Hcanop,   
     &   fwsoil,topfws,wscontent,           
     &   LAI,Sps,
     &   doy,hour,radsol,tair,dair,eairP,
     &   windU0,rain,wethr,
     &   Rnet,G,Esoil,Hcrop,Ecstot,Anet,
     &   Tsoil,DepH2O,
     &   wsmax,wsmin,  
     &   lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,
     &   stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,
     &   Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,
     &   H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,
     &   conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,
     &   Edjm,Entrpy,gam0,gam1,gam2,TminV,TmaxV,ToptV)

       real lat
       real gpp,evap,transp,LAI
       real tauL(3),rhoL(3),rhoS(3),reffbm(3),reffdf(3)
       real extkbm(3),extkdm(3)
       real Radabv(2),Qcan(3,2),Qcan0(3)
       real Acan(2),Ecan(2),Hcan(2),Tcan(2), Gbwcan(2), Gswcan(2)

       real topfws        
       integer idoy,ihour,ileaf
       integer jrain,i,j,k


       real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5),
     &   GbwcL(5),GswcL(5),hG(5),hIL(5)
       real Gaussx(5),Gaussw(5),Gaussw_cum(5)
      
       character*80 commts



       data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
       data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
       data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/


       call  yrday(doy,hour,lat,radsol,fbeam)
                

       jrain=int(wethr)
       if(jrain.lt.0) goto 19
       idoy=int(doy)
       hours=idoy*1.0+hour/24.0
       coszen=sinbet(doy,lat,pi,hour)             

       if(windU0.lt.0.01) windU0=0.01

       if(topfws.gt.0.5) then
         rhoS(2)=0.18
       else
         rhoS(2)=0.52-0.68*topfws
       endif


       FLAIT =LAI 
       radabv(1)=0.5*radsol                 
       radabv(2)=0.5*radsol                 

       Acanop=0.0
       Ecanop=0.0
       Hcanop=0.0
       fslt=0.0
       fsltx=0.0



       call xlayers(Sps,Tair,Dair,radabv,G,Esoil,fbeam,eairP,
     &   windU0,co2ca,fwsoil,LAI,coszen,idoy,hours,
     &   tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,
     &   Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,
     &   cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,
     &   gsw0,alpha,stom_n,
     &   Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &   extKb,
     &   Rnst1,Qcan1,Acan1,Ecan1,Hcan1,Gbwc1,Gswc1,Tleaf1,
     &   Rnst2,Qcan2,Acan2,Ecan2,Hcan2,Gbwc2,Gswc2,Tleaf2,
     &   Rcan1,Rcan2,Rsoilabs,Hsoil,
     &   RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,
     &   TminV,TmaxV,ToptV)


       do ng=1,5  
         hG(ng)=Gaussw_cum(ng)                   
         hIL(ng)=0.2*ng                          
       enddo
           



       AcanL(1)=AcanL(1)+(hIL(1)-hG(1))*(AcanL(2)-AcanL(1))
     &   /(hG(2)-hG(1))
       AcanL(2)=AcanL(2)+(hIL(2)-hG(2))*(AcanL(3)-AcanL(2))
     &   /(hG(3)-hG(2))
       AcanL(3)=AcanL(2)+(hIL(3)-hG(2))*(AcanL(3)-AcanL(2))
     &   /(hG(3)-hG(2))
       AcanL(4)=AcanL(3)+(hIL(4)-hG(3))*(AcanL(4)-AcanL(3))
     &   /(hG(4)-hG(3))
       EcanL(1)=EcanL(1)+(hIL(1)-hG(1))*(EcanL(2)-EcanL(1))
     &   /(hG(2)-hG(1))
       EcanL(2)=EcanL(2)+(hIL(2)-hG(2))*(EcanL(3)-EcanL(2))
     &   /(hG(3)-hG(2))
       EcanL(3)=EcanL(2)+(hIL(3)-hG(2))*(EcanL(3)-EcanL(2))
     &   /(hG(3)-hG(2))
       EcanL(4)=EcanL(3)+(hIL(4)-hG(3))*(EcanL(4)-EcanL(3))
     &   /(hG(4)-hG(3))
       HcanL(1)=HcanL(1)+(hIL(1)-hG(1))*(HcanL(2)-HcanL(1))
     &   /(hG(2)-hG(1))
       HcanL(2)=HcanL(2)+(hIL(2)-hG(2))*(HcanL(3)-HcanL(2))
     &   /(hG(3)-hG(2))
       HcanL(3)=HcanL(2)+(hIL(3)-hG(2))*(HcanL(3)-HcanL(2))
     &   /(hG(3)-hG(2))
       HcanL(4)=HcanL(3)+(hIL(4)-hG(3))*(HcanL(4)-HcanL(3))
     &   /(hG(4)-hG(3))



         mol_mass=44./1000.
         do ng=1,5
            EcanL(ng)=EcanL(ng)+Esoil
            HcanL(ng)=HcanL(ng)+Hsoil
          end do


         fsltx=(1.0-exp(-LAI*extkb))/(extkb*LAI)  


       Acan1=Acan1*1.0e6                     
       Acan2=Acan2*1.0e6
       Acanop=(Acan1+Acan2)       
       Qcan1=Qcan1                           
       Qcan2=Qcan2                 
       Qcanop=Qcan1+Qcan2 
       Ecanop=Ecan1+Ecan2                    
       Tcanop=Tleaf1*fsltx+Tleaf2*(1.0-fsltx)
       Hcanop=Hcan1+Hcan2                    
       if(Rsoilabs.LT.0.0) Rsoilabs=0.0
       if(Ecanop.LT.0.0)Ecanop=0.0

       gpp=Acanop*3600.0*12.0/1.0e6          
       transp=Ecanop*3600.0/((2.501-0.00236*Tair)*1000000.0)  
       if(transp.lt.0.0)transp=0.0
       evap=0.9*Rsoilabs*3600.0/2260.0/1000.0     
       if(evap.lt.0.0)evap=0.0

19     continue
       return
       end
c============================================================================

       subroutine respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,
     &   LAI,SLA,bmstem,bmroot,bmleaf,
     &   StemSap,RootSap,NSC,fnsc,
     &   RaLeaf,RaStem,RaRoot,Rauto)


       implicit none
       real LAIMIN,LAI,GPP,SLA
       real Tair,Tsoil,DepH2O
       real bmstem,bmroot,bmleaf,StemSap,RootSap
       real NSC,fnsc
       real Q10
       real RaLeaf,RaStem,RaRoot,Rauto
       real Rl0,Rs0,Rr0
       real c                  

       c=3600.*12./1000000.    
       Q10=2.0
       if(LAI.gt.LAIMIN) then
         Rl0=3.              
         Rs0=1.5
         Rr0=3.
         RaLeaf=Rl0*bmleaf*0.45*SLA*0.1*Q10**((Tair-25.)/10.)*fnsc*c
         RaStem=Rs0*StemSap*0.001 *Q10**((Tair-25.)/10.)*fnsc*c
         RaRoot=Rr0*RootSap*0.001 *Q10**((Tair-25.)/10.)*fnsc*c
       else
         RaLeaf=0.3*GPP
         RaStem=0.3*GPP
         RaRoot=0.4*GPP
       endif
       Rauto=Raleaf+Rastem+Raroot
       if(Rauto > 0.1*NSC)then
         Raleaf=Raleaf/Rauto*0.1*NSC
         Rastem=Rastem/Rauto*0.1*NSC
         Raroot=Rastem/Rauto*0.1*NSC
         Rauto=0.1*NSC
       endif
       return
       end


       subroutine soilwater(wsmax,wsmin,rdepth,FRLEN,
     &   rain,tair,transp,wcl,tsoil,Rh,thksl,LAI,       
     &   evap,runoff,wscontent,fwsoil,topfws,  
     &   omega,omega_S,WaterR,WaterS,SatFracL)  


       implicit none

       real wsmax,wsmin,wsmaxL(10),wsminL(10) 
       real FLDCAP,WILTPT,FLDCAPL(10),WILTPTL(10) 

       real LAI,rdepth
       integer nfr

       real precp,rain 
       real tair,TSOIL,ts          

       real evap,transp,evaptr,TEVAP,AEVAP

       real wscontent,fwsoil,topfws,omega,topomega,omega_S
       real fw(10),ome(10),W_signal
       real WaterR,WaterS(10),SatFracL(10)

       real RAWCL(10) 
       real thksl(10),depth(10),wsc(10),WUPL(10),EVAPL(10),SRDT(10)
       real plantup(10)
       real Tsrdt
       real frlen(10) 
       real wcl(10) 
       real fwcln(10) 
       real wtdeficit(10),DWCL(10),Tr_ratio(10)
       real Twater,Twater1,Twater2,Tthk,dWaterS,netflux
       real wtneed,wtadd,twtadd,infilt,runoff,roff_layer,tr_allo

       real RH,Rsoil,Rd,density,sp_heat,psychro,la,P
       real esat
       real exchangeL,supply,demand,omegaL(10)
       integer i,j,k

       WILTPT =wsmin/100.0
       FLDCAP =wsmax/100.0
       WILTPTL=wsmin/100.0
       FLDCAPL=wsmax/100.0

       twater=0.0
       twater1=0.0
       twater2=0.0
       nfr=0

       precp=rain

       do i=1,10
         wtdeficit(i)=0.0
         dwcl(i)=0.0
         evapl(i)=0.0
         WUPL(i)=0.0
         SRDT(i)=0.0
         DEPTH(i)=0.0
       enddo



       DEPTH(1)=10.0
       DO i=2,10
         DEPTH(i)=DEPTH(i-1)+THKSL(i)
       enddo
       do i=1,10
         IF(rdepth.GT.DEPTH(i)) nfr=i+1
       enddo
       IF (nfr.GT.10) nfr=10
       do i=1,10
         if(FLDCAPL(i).gt.wcl(i))wtdeficit(i)=FLDCAPL(i)-wcl(i)
       enddo


       infilt=precp  

       TWTADD=0
       roff_layer=0.0
       do i=1,10
         IF(infilt.GT.0.0)THEN

           WTADD=AMIN1(INFILT,wtdeficit(i)*thksl(i)*10.0) 

           WCL(i)=(WCL(i)*(thksl(i)*10.0)+WTADD)/(thksl(i)*10.0)
           FWCLN(I)=WCL(I)       
           TWTADD=TWTADD+WTADD       
           INFILT=INFILT-WTADD 
         END IF

         if(infilt.GT.0.0)THEN
           roff_layer=roff_layer + INFILT*0.05*(i-1)
           INFILT=INFILT -  INFILT*0.05*(i-1)
         endif
       enddo

       if(precp.gt.0.0.and.wcl(1).gt.wcl(2))then
         supply=(wcl(1)-wcl(2))/3.0
         wcl(1)=wcl(1)-2.0*supply
         wcl(2)=wcl(2)+supply
       endif


       runoff=INFILT + roff_layer   


       do i=1,10
         wsc(i)=Amax1(0.00,(wcl(i)-wiltpt)*THKSL(i)*10.0)
         omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAPL(i)-WILTPT))
       enddo
       supply=0.0
       demand=0.0
       do i=1,9
         if(omegaL(i).gt.0.3)then
           supply=wsc(i)/360.0*omegaL(i)
           demand=(FLDCAPL(i)-wcl(i+1))*THKSL(i+1)*10.0/360.0
     &       *(1.0-omegaL(i+1))
           exchangeL=AMIN1(supply,demand)
           wsc(i)=wsc(i)- exchangeL
           wsc(i+1)=wsc(i+1)+ exchangeL
           wcl(i)=wsc(i)/(THKSL(i)*10.0)+wiltpt
           wcl(i+1)=wsc(i+1)/(THKSL(i+1)*10.0)+wiltpt
         endif
       enddo



       if(wcl(1).LT.wiltpt)then
         evap=0.0
       else
         Rsoil=10.1*exp(1.0/wcl(1))
         Rd=20.5 
         P=101325.0  
         density=1.204 
         la=(2.501-0.00236*Tair)*1000000.0 
         sp_heat=1012.0  
         psychro=1628.6*P/la

         evap=1.0*esat(tair)*(1.0-RH/100.0)/
     &     (Rsoil+Rd)*density*sp_heat/psychro/la*3600.0
       endif



       Twater=0
       do i=1,10
         wsc(i)=(wcl(i)-wiltpt)*THKSL(I)*10.0
         Twater=Twater+wsc(i)  
       enddo

       Tsrdt=0.0
       do i=1,10

         SRDT(I)=EXP(-4.73*(DEPTH(I)-THKSL(I)/2.0)/100.0)/1.987
         Tsrdt=Tsrdt+SRDT(i)/(i*i)  
       enddo

       do i=1,10
         SRDT(i)=SRDT(i)/Tsrdt
       enddo

       do i=1,10
         EVAPL(I)=Amax1(AMIN1(evap*SRDT(i),wsc(i)),0.0)  
         DWCL(I)=EVAPL(I)/(THKSL(I)*10.0) 
       enddo


       do i=1,10
         wcl(i)=wcl(i)-DWCL(i)
       enddo

       evap=0.0       
       do i=1,10
         evap=evap+EVAPL(I)
       enddo

       Twater=0
       do i=1,nfr
         wsc(i)=(wcl(i)-wiltpt)*THKSL(I)*10.0
         Twater=Twater+AMAX1(wsc(i),0.0) 
       enddo
       if(transp.gt.Twater/2.0)transp=Twater/2.0                     
       tr_allo=0.0
       do i=1,nfr
         tr_ratio(i)=FRLEN(i) 
         tr_allo=tr_allo+tr_ratio(i)
       enddo

       do i=1,nfr
         plantup(i)=AMIN1(transp* tr_ratio(i)/tr_allo, wsc(i)) 
         wupl(i)=plantup(i)/(thksl(i)*10.0)
         wcl(i)=wcl(i)-wupl(i)
       end do

       transp=0.0
       do i=1,nfr
         transp=transp+plantup(i)
       enddo


       Twater=0
       Tthk=0
       do i=1,nfr
         Twater=Twater+wcl(i)*THKSL(I)*10.0 
         Tthk=Tthk+thksl(i)*10.0 
       enddo
       wscontent=Twater/Tthk
       if(wscontent.lt.WILTPT) wscontent=WILTPT+0.00001
       omega_S=(wscontent-WILTPT)/(FLDCAP-WILTPT)
       fwsoil=amin1(1.0,3.333*omega)
       topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))
       if(fwsoil.lt.0.0) fwsoil=0.000001
       if(topfws.lt.0.0) topfws=0.000001
       if(omega.lt.0.0) omega=0.0000001
       Twater=Twater-WILTPT*Tthk
       WaterR=Twater+WILTPT*Tthk
      

       do i=1,10 
         ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
         WaterS(i)=wcl(i)*THKSL(I)*10.0
         ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
         SatFracL(i)=ome(i)
         fw(i)=amin1(1.0,3.333*ome(i))
       enddo
       fwsoil=0.0
       omega=0.0
       do i=1,nfr
         fwsoil=fwsoil+fw(i)*frlen(i)
         omega=omega+ome(i)*frlen(i)
       enddo
       return
       end




       subroutine plantgrowth(Tair,Tavg72,omega,GLmax,GRmax,GSmax,
     &   LAI,LAIMAX,LAIMIN,SLA,Tau_L,
     &   bmleaf,bmroot,bmstem,bmplant,
     &   Rootmax,Stemmax,SapS,
     &   StemSap,RootSap,Storage,GDD5,
     &   stor_use,onset,accumulation,gddonset,
     &   Sps,NSC,fnsc,NSCmin,NSCmax,
     &   store,add,L_fall,Tcold,Gamma_Wmax,Gamma_Tmax,
     &   NPP,alpha_L,alpha_W,alpha_R)
       implicit none
       real NSC,NSCmin,NSCmax,fnsc
       real store,Storage,GDD5,stor_use,accumulation,gddonset
       integer onset,duration,offset,dormancy
       real GLmax,GRmax,GSmax,TauLeaf
       real GrowthP,GrowthL,GrowthR,GrowthS
       real Tair,Tavg72,T72(72)
       real omega,LAI,LAIMAX,LAIMIN,SLA

       real bmleaf,bmroot,bmstem,bmplant,NPP
       real Rootmax,Stemmax,SapS,SapR
       real bmL,bmR,bmP,bmS,StemSap,RootSap

       real St,Sw,Ss,Sn,SL_rs,SR_rs,Slai,Sps
       real RS,RS0,RSw
       real gamma_W,gamma_Wmax,gamma_T,gamma_Tmax,gamma_N
       real beta_T,Tcold,Twarm
       real bW,bT,W
       real L_fall,L_add,add,NL_fall,NL_add,Tau_L
       real alpha_L,alpha_W,alpha_R,alpha_St
       real Twsoil(7),Tavg
       integer i

       bmL=bmleaf
       bmR=bmRoot
       bmS=bmStem
       bmP=bmPlant

       bW=2.0
       bT=2.0

       Twarm=30.0

       RS0=1.0
       RS=bmR/bmL

       if(bmL.lt.NSC/0.333*0.5)bmL=NSC/0.333*0.5
       if(bmR.lt.NSC/0.333*0.5)bmR=NSC/0.333*0.5
       if(bmS.lt.NSC/0.334*0.5)bmS=NSC/0.334*0.5
       StemSap=MIN(Stemmax,SapS*bmS)  
       RootSap=bmR
       if(StemSap.lt.0.001)StemSap=0.001
       if(RootSap.lt.0.001)RootSap=0.001
       beta_T=1


       if((GDD5.gt.gddonset).and.onset.eq.0.and.storage.gt.stor_use)then
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
         store=0.005*NSC
       else
         store=0.0
       endif
       accumulation=accumulation+store

       Sps=1.0 - fnsc
       sps=AMAX1(0.001,sps)
       St=1./(1.+19.*EXP(-0.2*(Tair)))
       Sw=AMAX1(0.333, 0.333+omega)
       W=AMIN1(1.0,50.*omega)
       Ss=AMIN1(1.0,2.*fnsc)

       SL_rs=RS/(RS+RS0*(1.5-omega))
       SR_rs=(RS0*(1.5-omega))/(RS+RS0*(1.5-omega))
       Slai=amin1(1.0,2.333*(LAIMAX-LAI)/(LAIMAX-LAIMIN))

       GrowthL=GLmax*bmL    *St*Sw*fnsc*SL_rs*Slai*0.45
       GrowthR=GRmax*RootSap*St*Sw*fnsc*SR_rs*0.5 *0.45       
       GrowthS=GSmax*StemSap*St*Sw*fnsc*0.5*0.5 *0.45

       if(GrowthL.LT.0.0)GrowthL=0.0
       if(GrowthR.LT.0.0)GrowthR=0.0
       if(GrowthS.LT.0.0)GrowthS=0.0

       GrowthP=GrowthL + GrowthR + GrowthS
       if(GrowthP.gt.NSC*0.5)then
         GrowthL=0.5*NSC*GrowthL/GrowthP
         GrowthR=0.5*NSC*GrowthR/GrowthP
         GrowthS=0.5*NSC*GrowthS/GrowthP
       endif
       NPP = add + GrowthL + GrowthR + GrowthS
       if(NPP.eq.0.0)then
         alpha_L=0.333
         alpha_W=0.333
         alpha_R=0.333
       else
         alpha_L=(GrowthL+add)/NPP
         alpha_W=GrowthS/NPP
         alpha_R=GrowthR/NPP
       endif

       if(Tair.gt.(Tcold+10.)) then
         beta_T=1.
       else 
         if(Tair.gt.Tcold)beta_T=(Tair-Tcold)/10.
         if(Tair.LE.Tcold)beta_T=0.
       endif
       if (tau_L < 8760)then
         gamma_W=(1. - W)     **bW * gamma_Wmax
         gamma_T=(1. - beta_T)**bT * gamma_Tmax
         gamma_N=0.0 
       else
         gamma_W=0.
         gamma_T=0.
         gamma_N=1.0/(Tau_L*Sw)
       endif

       if(LAI < LAIMIN) then
         gamma_W=0.
         gamma_T=0.
         gamma_N=0.
       endif
       L_fall=bmleaf*0.45*AMIN1((gamma_W+gamma_T+gamma_N),0.99)
       return
       end



       subroutine TCS(Tair,Tsoil,omega,
     &   NPP,alpha_L,alpha_W,alpha_R,
     &   L_fall,tau_L,tau_W,tau_R,
     &   tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass,
     &   Q_leaf,Q_wood,Q_root,
     &   Q_fine,Q_coarse,Q_Micr,Q_Slow,Q_Pass,
     &   Rh_f,Rh_c,Rh_Micr,Rh_Slow,Rh_Pass,S_w_min,Q10_h)
       implicit none
       real NPP,NPP_L,NPP_W,NPP_R
       real L_fall,L_add,LAI,SLA
       real Tair,Tsoil,omega

       real alpha_L,alpha_W,alpha_R

       real tau_L,tau_W,tau_R
       real tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass

       real Q_leaf,Q_wood,Q_root
       real Q_fine,Q_coarse,Q_Micr,Q_Slow,Q_Pass
       real eta 
       real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M

       real f_CO2_fine,f_CO2_coarse,f_CO2_Micr,f_CO2_Slow,f_CO2_Pass

       real tau_L_a,tau_W_a,tau_R_a
       real tau_F_a,tau_C_a,tau_Micr_a,tau_Slow_a,tau_Pass_a

       real Out_leaf,Out_wood,Out_root
       real Out_fine,Out_coarse,Out_Micr,Out_Slow,Out_Pass

       real Rh_f,Rh_c,Rh_Micr,Rh_Slow,Rh_Pass,Q10_h

       real TminV,TmaxV,ToptV,Tcold,Gamma_Wmax,Gamma_Tmax

       real S_omega 
       real S_t     
       real Ta,Ts   
       real S_w_min    
       real Tref,T0,T,E0,mid
       integer i,j,k,n,m
       integer day,week,month,year




       S_omega=S_w_min + (1.-S_w_min)*amin1(1.0,2.0*omega)
       S_t=1.5-1./(1.+19.*exp(-0.15*(Tsoil-15.0)))


       NPP_L=alpha_L*NPP            
       NPP_W=alpha_W*NPP
       NPP_R=alpha_R*NPP





       S_t=Q10_h**((Tsoil-25.)*0.1)
       Out_leaf=L_fall
       Out_wood=Q_wood/tau_W    * S_T*S_omega
       Out_root=Q_root/tau_R    * S_T*S_omega
       Out_fine=Q_fine/tau_F    * S_T*S_omega
       Out_coarse=Q_coarse/tau_C* S_T*S_omega
       Out_Micr=Q_micr/tau_Micr * S_T*S_omega
       Out_Slow=Q_Slow/tau_Slow * S_T*S_omega
       Out_Pass=Q_Pass/tau_Pas s* S_T*S_omega


       eta=0.15 
       f_F2M=0.45
       f_C2M=0.275
       f_C2S=0.275
       f_M2S=0.296
       f_M2P=0.002
       f_S2P=0.015
       f_S2M=0.42
       f_P2M=0.45

       Q_leaf  =Q_leaf - Out_leaf + NPP_L    
       Q_wood  =Q_wood - Out_wood + NPP_W
       Q_root  =Q_root - Out_root + NPP_R

       Q_fine  =Q_fine - Out_fine + Out_leaf + eta*Out_wood + Out_root
       Q_coarse=Q_coarse - Out_coarse + (1.-eta)*Out_wood
       Q_Micr = Q_Micr - Out_Micr
     &   + f_F2M*Out_fine+f_C2M*Out_coarse
     &   + f_S2M*Out_Slow+f_P2M*Out_Pass
       Q_Slow = Q_Slow - Out_Slow
     &   + f_C2S*Out_coarse + f_M2S*Out_Micr
       Q_Pass = Q_Pass - Out_Pass
     &   + f_M2P*Out_Micr + f_S2P*Out_Slow

       Rh_f =  Out_fine   * (1. - f_F2M)
       Rh_c =  Out_coarse * (1. - f_C2M - f_C2S)
       Rh_Micr=Out_Micr   * (1. - f_M2S - f_M2P)
       Rh_Slow=Out_Slow   * (1. - f_S2P - f_S2M)
       Rh_Pass=Out_Pass   * (1. - f_P2M)

       return
       end


       subroutine xlayers(Sps,Tair,Dair,radabv,G,Esoil,fbeam,eairP,
     &   windU0,co2ca,fwsoil,FLAIT,coszen,idoy,hours,
     &   tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,
     &   Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,
     &   cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,
     &   gsw0,alpha,stom_n,
     &   Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &   extKb,
     &   Rnst1,Qcan1,Acan1,Ecan1,Hcan1,Gbwc1,Gswc1,Tleaf1,
     &   Rnst2,Qcan2,Acan2,Ecan2,Hcan2,Gbwc2,Gswc2,Tleaf2,
     &   Rcan1,Rcan2,Rsoilabs,Hsoil,
     &   RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,
     &   TminV,TmaxV,ToptV)









       real Gaussx(5),Gaussw(5)
       real layer1(5),layer2(5)
       real tauL(3),rhoL(3),rhoS(3),Qabs(3,2),Radabv(2),Rnstar(2)
       real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
       real gbleaf(2),gsleaf(2),QSabs(3,2),Qasoil(2)
       integer ng,nw
       real rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       


       real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5),
     &   GbwcL(5),GswcL(5)
      


       data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
       data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/


       Rnst1=0.0        
       Rnst2=0.0        
       Qcan1=0.0        
       Qcan2=0.0
       Rcan1=0.0        
       Rcan2=0.0
       Acan1=0.0        
       Acan2=0.0
       Ecan1=0.0        
       Ecan2=0.0
       Hcan1=0.0        
       Hcan2=0.0
       Gbwc1=0.0        
       Gbwc2=0.0
       Gswc1=0.0        
       Gswc2=0.0
       Tleaf1=0.0       
       Tleaf2=0.0  
  

       raero=50./windU0                           


       xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
       xphi2 = 0.877 * (1.0 - 2.0*xphi1)
       funG=xphi1 + xphi2*coszen                             
      
       if(coszen.gt.0) then                                  
         extKb=funG/coszen                                   
       else
         extKb=100.
       end if




       pi180=3.1416/180.
       cozen15=cos(pi180*15)
       cozen45=cos(pi180*45)
       cozen75=cos(pi180*75)
       xK15=xphi1/cozen15+xphi2
       xK45=xphi1/cozen45+xphi2
       xK75=xphi1/cozen75+xphi2
       transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+
     &   0.178*exp(-xK75*FLAIT)
       extkd=(-1./FLAIT)*alog(transd)
       extkn=extkd             





       do nw=1,2              
       scatt(nw)=tauL(nw)+rhoL(nw)                      
       if((1.-scatt(nw))<0.0)scatt(nw)=0.9999           
       kpr(nw,1)=extKb*sqrt(1.-scatt(nw))               
       kpr(nw,2)=extkd*sqrt(1.-scatt(nw))             
       rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            
       rhoc15=2.*xK15*rhoch/(xK15+extkd)                                
       rhoc45=2.*xK45*rhoch/(xK45+extkd)
       rhoc75=2.*xK75*rhoch/(xK75+extkd)
       rhoc(nw,2)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
       rhoc(nw,1)=2.*extKb/(extKb+extkd)*rhoch                          
       reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))                      
     &   *exp(-2.*kpr(nw,1)*FLAIT) 
       reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))                      
     &   *exp(-2.*kpr(nw,2)*FLAIT)  
       enddo



       call Radiso(flait,flait,Qabs,extkd,Tair,eairP,cpair,Patm,
     &   fbeam,airMa,Rconst,sigma,emleaf,emsoil,
     &   emair,Rnstar,grdn)

       TairK=Tair+273.2


       do ng=1,5
         flai=gaussx(ng)*FLAIT

         call goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,
     &     scatt,xfang,Qabs) 

         call Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,
     &     fbeam,airMa,Rconst,sigma,emleaf,emsoil,
     &     emair,Rnstar,grdn)
         windUx=windU0*exp(-extkU*flai)             
         scalex=exp(-extkn*flai)                    
         Vcmxx=Vcmx0*scalex
         eJmxx=eJmx0*scalex
         if(radabv(1).ge.10.0) then                          

           call agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,
     &       co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,
     &       Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,
     &       gsw0,alpha,stom_n,
     &       Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,
     &       Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &       Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,
     &       TminV,TmaxV,ToptV)
         else
           call agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,
     &       co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,
     &       Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,
     &       gsw0,alpha,stom_n,
     &       Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,
     &       Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &       Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
         endif  
         fslt=exp(-extKb*flai)                        
         fshd=1.0-fslt                                
         Rnst1=Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  
         Rnst2=Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
         RnstL(ng)=Rnst1+Rnst2

         Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  
         Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
         QcanL(ng)=Qcan1+Qcan2

         Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  
         Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
         RcanL(ng)=Rcan1+Rcan2

         if(Aleaf(1).lt.0.0)Aleaf(1)=0.0      
         if(Aleaf(2).lt.0.0)Aleaf(2)=0.0      

         Acan1=Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    
         Acan2=Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
         AcanL(ng)=Acan1+Acan2

         layer1(ng)=Aleaf(1)
         layer2(ng)=Aleaf(2)

         Ecan1=Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
         Ecan2=Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
         EcanL(ng)=Ecan1+Ecan2

         Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
         Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
         HcanL(ng)=Hcan1+Hcan2

         Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n

         Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n

         Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
         Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT

200      continue
       enddo

       FLAIT1=(1.0-exp(-extKb*FLAIT))/extkb
       Tleaf1=Tleaf1/FLAIT1
       Tleaf2=Tleaf2/(FLAIT-FLAIT1)


       Rsoilab1=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)
     &   +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          
       Rsoilab2=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)
     &   +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          
       Rsoilab1=Rsoilab1*Radabv(1)
       Rsoilab2=Rsoilab2*Radabv(2)

       Tlk1=Tleaf1+273.2
       Tlk2=Tleaf2+273.2
       QLair=emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
       QLleaf=emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)
     &   +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
       QLleaf=QLleaf*(1.0-exp(-extkd*FLAIT)) 
       QLsoil=emsoil*sigma*(TairK**4)
       Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil






       Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3 



       Esoil=0.9*(Rsoilabs-G)
       Hsoil=0.1*(Rsoilabs-G)

       return
       end 


       subroutine goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,
     &   scatt,xfang,Qabs)
     

















       real radabv(2)
       real Qabs(3,2),reff(3,2),kpr(3,2),scatt(2)
       xu=coszen                                         
      

       xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
       xphi2 = 0.877 * (1.0 - 2.0*xphi1)
       funG=xphi1 + xphi2*xu                             
      
       if(coszen.gt.0) then                                  
         extKb=funG/coszen                                   
       else
         extKb=100.
       end if
                       

       do nw=1,2
         Qd0=(1.-fbeam)*radabv(nw)                                          
         Qb0=fbeam*radabv(nw)
         Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*
     &     FLAI))+Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-
     &     extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
         Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))
       end do
       return
       end


       subroutine Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,
     &   fbeam,airMa,Rconst,sigma,emleaf,emsoil,
     &   emair,Rnstar,grdn)




       real Rnstar(2)
       real Qabs(3,2)
       TairK=Tair+273.2


       rhocp=cpair*Patm*airMa/(Rconst*TairK)   


       emsky=0.642*(eairP/Tairk)**(1./7)       
     

       ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)
       tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            
       emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      




       emair=emsky

       if(emair.gt.1.0) emair=1.0
      


       Bn0=sigma*(TairK**4.)
       Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)
     &   + exp(-extkd*(flait-flai))*(emsoil-emleaf))


       Rnstar(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
       Rnstar(2)=Qabs(1,2)+Qabs(2,2)+Bnxi

       grdn=4.*sigma*(TairK**3.)*extkd*emleaf*
     &   *(exp(-extkd*flai)+exp(-extkd*(flait-flai)))
     &   /rhocp
       return
       end

       subroutine agsean_day(Sps,Qabs,Rnstar,grdn,windUx,Tair,
     &   Dair,co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,
     &   hours,Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,
     &   H2OMw,Dheat,gsw0,alpha,stom_n,
     &   Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &   Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci,
     &   TminV,TmaxV,ToptV)

       integer kr1,ileaf
       real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
       real gbleaf(2), gsleaf(2)
       real Qabs(3,2),Rnstar(2)

       TairK=Tair+273.2
       rhocp=cpair*Patm*AirMa/(Rconst*TairK)
       H2OLv=H2oLv0-2.365e3*Tair
       slope=(esat(Tair+0.1)-esat(Tair))/0.1
       psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
       Cmolar=Patm/(Rconst*TairK)
       weighJ=1.0


       if(windUx/wleaf>=0.0)then
         gbHu=0.003*sqrt(windUx/wleaf)    
       else
         gbHu=0.003 
       endif         
       do ileaf=1,2              

         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    

         Dleaf=Dair                

         co2cs=co2ca               
         Qapar = (4.6e-6)*Qabs(1,ileaf)

         kr1=0                     

         do               

           Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3.)     
           gbHf=0.5*Dheat*(Gras**0.25)/wleaf
           gbH=gbHu+gbHf                         
           rbH=1./gbH                            
           rbw=0.93*rbH                          

           rbH_L=rbH*stom_n/2.                   
           rrdn=1./grdn
           Y=1./(1.+ (rbH_L+raero)/rrdn)

           gbc=Cmolar*gbH/1.32            
           gsc0=gsw0/1.57                 
           varQc=0.0
           weighR=1.0
           call photosyn(Sps,CO2Ca,CO2Csx,Dleaf,Tlk,Qapar,Gbc, 
     &       theta,a1,Ds0,fwsoil,varQc,weighR,
     &       gsc0,alpha,Vcmxx,eJmxx,weighJ,
     &       conKc0,conKo0,Ekc,Eko,o2ci,Rconst,Trefk,
     &       Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &       Aleafx,Gscx,TminV,TmaxV,ToptV)  

           Aleaf(ileaf) = Aleafx      

           co2cs = co2ca-Aleaf(ileaf)/gbc
           co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc0

           gsw=gsc0*1.56       
           gswv=gsw/Cmolar                           
           rswv=1./gswv

           Eleaf(ileaf)=9.0*
     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    
     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))

           Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))

           Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp

           Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
           gbleaf(ileaf)=gbc*1.32*1.075
           gsleaf(ileaf)=gsw

           if(abs(Tlk1-Tlk).le.0.1) exit 

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
           endif                     
         enddo
       enddo
       return
       end

       subroutine agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,
     &   co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,
     &   Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,
     &   gsw0,alpha,stom_n,
     &   Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &   Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
       integer kr1,ileaf
       real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
       real gbleaf(2), gsleaf(2)
       real Qabs(3,2),Rnstar(2)

       TairK=Tair+273.2
       rhocp=cpair*Patm*AirMa/(Rconst*TairK)
       H2OLv=H2oLv0-2.365e3*Tair
       slope=(esat(Tair+0.1)-esat(Tair))/0.1
       psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
       Cmolar=Patm/(Rconst*TairK)
       weighJ=1.0



       gbHu=0.003*sqrt(windUx/wleaf)    

       do ileaf=1,2                  

         Tleaf(ileaf)=Tair
         Tlk=Tleaf(ileaf)+273.2    

         Dleaf=Dair                

         co2cs=co2ca               
         Qapar = (4.6e-6)*Qabs(1,ileaf)

         kr1=0                     
         do

           Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     
           gbHf=0.5*Dheat*(Gras**0.25)/wleaf
           gbH=gbHu+gbHf                         
           rbH=1./gbH                            
           rbw=0.93*rbH                          

           rbH_L=rbH*stom_n/2.                   
           rrdn=1./grdn
           Y=1./(1.+ (rbH_L+raero)/rrdn)

           gbc=Cmolar*gbH/1.32            
           gsc0=gsw0/1.57                        
           varQc=0.0                  
           weighR=1.0

           Aleafx=-0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
           gsc=gsc0

           Aleaf(ileaf) = Aleafx                     

           co2cs = co2ca-Aleaf(ileaf)/gbc
           co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc

           gsw=gsc*1.56                              
           gswv=gsw/Cmolar                           
           rswv=1./gswv

           Eleaf(ileaf)=
     &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/
     &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))

           Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))

           Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp

           Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
           gbleaf(ileaf)=gbc*1.32*1.075
           gsleaf(ileaf)=gsw


           if(abs(Tlk1-Tlk).le.0.1)exit
           if(kr1.gt.1000)exit

           Tlk=Tlk1 
           Tleaf(ileaf)=Tlk1-273.2
           kr1=kr1+1
         enddo                          
10       continue
       enddo
       return
       end

       subroutine ciandA(Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad)

       b2=g0+X*(Gma-Rd)
       b1=(1.-co2cs*X)*(Gma-Rd)+g0*(Bta-co2cs)-X*(Gma*gammas+Bta*Rd)
       b0=-(1.-co2cs*X)*(Gma*gammas+Bta*Rd)-g0*Bta*co2cs

       bx=b1*b1-4.*b2*b0
       if(bx.gt.0.0)then 

         ciquad = (-b1+sqrt(bx))/(2.*b2)
       endif

       if(ciquad.lt.0.or.bx.lt.0.)then
         Aquad = 0.0
         ciquad = co2Cs
       else
         Aquad = Gma*(ciquad-gammas)/(ciquad+Bta)
       end if
       return
       end


       subroutine goud1(FLAIT,coszen,radabv,fbeam,
     &   Tair,eairP,emair,emsoil,emleaf,sigma,
     &   tauL,rhoL,rhoS,xfang,extkb,extkd,
     &   reffbm,reffdf,extkbm,extkdm,Qcan)



















       integer nW
       real radabv(3)
       real rhocbm(3),rhocdf(3)
       real reffbm(3),reffdf(3),extkbm(3),extkdm(3)
       real tauL(3),rhoL(3),rhoS(3),scatL(3)
       real Qcan(3,2),Qcan0(3)


       fdiff=1.0-fbeam
       xu=coszen
       xphi1 = 0.5 -0.633*xfang - 0.33*xfang*xfang
       xphi2 = 0.877 * (1.0 - 2.0*xphi1)
       funG = xphi1 + xphi2*xu
       extkb=funG/xu
                       

       pi180=3.1416/180.
       cozen15=cos(pi180*15)
       cozen45=cos(pi180*45)
       cozen75=cos(pi180*75)
       xK15=xphi1/cozen15+xphi2
       xK45=xphi1/cozen45+xphi2
       xK75=xphi1/cozen75+xphi2
       transd=0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+
     &   0.178*exp(-xK75*FLAIT)
       extkd=(-1./FLAIT)*alog(transd)


       do nw=1,2
         scatL(nw)=tauL(nw)+rhoL(nw)
         if((1.-scatL(nw))<0.0) scatL(nw)=0.9999
         extkbm(nw)=extkb*sqrt(1.-scatL(nw))
         extkdm(nw)=extkd*sqrt(1.-scatL(nw))
         rhoch=(1.-sqrt(1.-scatL(nw)))/(1.+sqrt(1.-scatL(nw)))

         rhoc15=2.*xK15*rhoch/(xK15+extkd)
         rhoc45=2.*xK45*rhoch/(xK45+extkd)
         rhoc75=2.*xK75*rhoch/(xK75+extkd)
       
         rhocbm(nw)=2.*extkb/(extkb+extkd)*rhoch
         rhocdf(nw)=0.308*rhoc15+0.514*rhoc45+0.178*rhoc75

         reffbm(nw)=rhocbm(nw)+(rhoS(nw)-rhocbm(nw))
     &             *exp(-2.*extkbm(nw)*FLAIT)
         reffdf(nw)=rhocdf(nw)+(rhoS(nw)-rhocdf(nw))
     &             *exp(-2.*extkdm(nw)*FLAIT)  


         abshdn=fdiff*(1.0-reffdf(nw))*extkdm(nw)
     &     *(funE(extkdm(nw),FLAIT)-funE((extkb+extkdm(nw)),FLAIT))
     &     +fbeam*(1.0-reffbm(nw))*extkbm(nw)

     &     *(funE(extkbm(nw),FLAIT)-funE((extkb+extkbm(nw)),FLAIT))
     &     -fbeam*(1.0-scatL(nw))*extkb
     &     *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))

         absltn=fdiff*(1.0-reffdf(nw))*extkdm(nw)                             
     &     *funE((extkb+extkdm(nw)),FLAIT)                         
     &     +fbeam*(1.0-reffbm(nw))*extkbm(nw)

     &     *funE((extkb+extkbm(nw)),FLAIT)
     &     +fbeam*(1.0-scatL(nw))*extkb
     &     *(funE(extkb,FLAIT)-funE(2.0*extkb,FLAIT))



         Qcan(nw,1)=absltn*radabv(nw)

         Qcan(nw,2)=abshdn*radabv(nw)
       enddo


       TairK=Tair+273.2
      

       emsky=0.642*(eairP/Tairk)**(1./7)      


       ep8z=0.24+2.98e-12*eairP*eairP*exp(3000.0/TairK)
       tau8=amin1(1.0,1-ep8z*(1.4-0.4*ep8z))                
       emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4        




       emair=emsky
       if(emair.gt.1.0) emair=1.0                             

       Bn0=sigma*(TairK**4)
       QLW1=-extkd*emleaf*(1.0-emair)*funE((extkd+extkb),FLAIT)
     &   -extkd*(1.0-emsoil)*(emleaf-emair)*exp(-2.0*extkd*FLAIT)
     &   *funE((extkb-extkd),FLAIT)
       QLW2=-extkd*emleaf*(1.0-emair)*funE(extkd,FLAIT)
     &   -extkd*(1.0-emsoil)*(emleaf-emair)
     &   *(exp(-extkd*FLAIT)-exp(-2.0*extkd*FLAIT))/extkd
     &   -QLW1
       Qcan(3,1)=QLW1*Bn0
       Qcan(3,2)=QLW2*Bn0
       return
       end


       subroutine photosyn(Sps,CO2Ca,CO2Csx,Dleafx,Tlkx,Qaparx,Gbcx,
     &   theta,a1,Ds0,fwsoil,varQc,weighR,
     &   g0,alpha,
     &   Vcmx1,eJmx1,weighJ,conKc0,conKo0,Ekc,Eko,o2ci,
     &   Rconst,Trefk,Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,
     &   Aleafx,Gscx,TminV,TmaxV,ToptV)






       if(Qaparx.le.0.) then                            
         Aleafx=-0.0089*Vcmx1*exp(0.069*(Tlkx-293.2))   
         Gscx=g0
       endif

      
       TminJ=TminV
       TmaxJ=TmaxV
       ToptJ=ToptV 
      
       Tlf=Tlkx-273.2
       VcmxT=VJtemp(Tlf,TminV,TmaxV,ToptV,Vcmx1)
       eJmxT=VJtemp(Tlf,TminJ,TmaxJ,ToptJ,eJmx1)      

       eJ = weighJ*fJQres(eJmxT,alpha,Qaparx,theta)

       conKcT = EnzK(Tlkx,Trefk,conKc0,Rconst,Ekc)
       conKoT = EnzK(Tlkx,Trefk,conKo0,Rconst,Eko)


       Rd = 0.0089*VcmxT*weighR                              
       Tdiff=Tlkx-Trefk
       gammas = gam0*(1.+gam1*Tdiff+gam2*Tdiff*Tdiff)       

       gamma = 0.0




       X = a1*fwsoil/((co2csx - gamma)*(1.0 + Dleafx/Ds0))

       Gma = VcmxT  
       Bta = conKcT*(1.0+ o2ci/conKoT)
       call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci2,Acx)

       Gma = eJ/4.
       Bta = 2.*gammas

       call ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,co2ci4,Aqx)

       sps=AMAX1(0.001,sps)                  
       Aleafx = (amin1(Acx,Aqx) - Rd)*sps     

       CO2csx = co2ca-Aleafx/Gbcx
       Gscx=g0+X*Aleafx  

       return
       end

       function funeJ(alpha,eJmxT,Qaparx)
       funeJ=alpha*Qaparx*eJmxT/(alpha*Qaparx+2.1*eJmxT)
       return
       end

       real function esat(T)

       esat=610.78*exp(17.27*T/(T+237.3))
       return
       end


       real function evapor(Td,Tw,Patm)

       gamma = (64.6 + 0.0625*Td)/1.e5
       evapor = esat(Tw)- gamma*(Td-Tw)*Patm
       return
       end


       real function Vjmax(Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
       anum = Vjmax0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
       aden = 1. + EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
       Vjmax = anum/aden
       return
       end

       real function funE(extkbd,FLAIT)
       funE=(1.0-exp(-extkbd*FLAIT))/extkbd
       return
       end




       real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
       if(Tlf.lt.TminVJ) Tlf=TminVJ   
       if(Tlf.gt.TmaxVJ) Tlf=TmaxVJ
       pwr=(TmaxVJ-ToptVJ)/(ToptVj-TminVj)
       VJtemp=VJmax0*((Tlf-TminVJ)/(ToptVJ-TminVJ))*
     &   ((TmaxVJ-Tlf)/(TmaxVJ-ToptVJ))**pwr 
       return
       end


       real function fJQres(eJmx,alpha,Q,theta)
       AX = theta                                 
       BX = alpha*Q+eJmx                          
       CX = alpha*Q*eJmx                          
       if((BX*BX-4.*AX*CX)>=0.0)then
         fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
       else
         fJQres = (BX)/(2*AX)                   
       endif

       return
       end


       real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)

       temp1=(Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk)
       EnzK=EnzK0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
       return
       end


       real function sinbet(doy,lat,pi,timeh)
       real lat

       rad = pi/180.

       sinlat = sin(rad*lat)
       coslat = cos(rad*lat)

       sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
       cosdec=sqrt(1.-sindec*sindec)

       A = sinlat*sindec
       B = coslat*cosdec
       sinbet = A+B*cos(pi*(timeh-12.)/12.)
       return
       end


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
       if(tmprat.gt.0.22.and.tmprat.le.0.35)then
         fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
       endif
       if(tmprat.gt.0.35.and.tmprat.le.tmpK)then
         fdiff=1.47-1.66*tmprat
       endif
       if(tmprat.ge.tmpK) then
         fdiff=tmpR
       endif
       fbeam=1.0-fdiff
       if(fbeam.lt.0.0) fbeam=0.0
       return
       end

       SUBROUTINE init_random_seed()
         INTEGER :: i, n, clock
         INTEGER, DIMENSION(:), ALLOCATABLE :: seed
         
         CALL RANDOM_SEED(size = n)
         ALLOCATE(seed(n))
          
         CALL SYSTEM_CLOCK(COUNT=clock)
          
         seed = clock + 37 * (/ (i - 1, i = 1, n) /)
         CALL RANDOM_SEED(PUT = seed)
          
         DEALLOCATE(seed)
       END SUBROUTINE

