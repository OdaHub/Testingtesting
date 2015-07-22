!BM-----------------------------------------------------------
Module summary
contains

!BS------------------------------------------------------
subroutine make_summary
use stringlengths
implicit none
character(len=LSTR)::fsum

print *,"Make summary file"
fsum="summary.tex"
open(1,file="summary.tex")
write(1,*) "\begin{verbatim}"
  call write_results(fsum)
write(1,*) "\end{verbatim}"
close(1)
end subroutine make_summary 
!ES------------------------------------------------------

!BS-------------------------------------------------------------
subroutine write_results(fsum)
use var
use results
use erfc
implicit none
character(len=*)::fsum
double precision::c1,idetect,flux,re_flux,corflux,pMD1,pMD2,pMD3
double precision::corMD1,corMD2,corMD3,teffMD1,teffMD2,teffMD3,teffflux
double precision::totrelerrPcross,totsimPcross,totPcross
integer::i
double precision, allocatable::teffTIS(:),TISsim(:),TIScyclesOPT(:)
double precision::rate,relerrRate,totsimrate,teffpcross,teffrate
double precision::teffpcrossOPT,teffrateOPT
double precision::Ncycmin,cyclesOptmin
double precision::simT0,simT1,teffT0,teffT1,CyclesOptT0
double precision::fluxPPS,teffPPS,rePPS,cyclesoptflux,totsimflux
integer::fluxcycles
double precision::PxinATS,PxinATS2,PxinAL,PxinAM,PxinAR,kapTST,kappaTSTBrownian
double precision::kapTSTtheor,wb,tol,x

allocate(teffTIS(numint))
allocate(TISsim(numint))
allocate(TIScyclesOpt(numint))

write(1,*) "MD results:"
write(1,*) "flux1: flux through:",MD_INTFL
write(1,*) "flux2: effective flux through:",MD_INTFM
write(1,*) "flux3: effective flux through:",MD_INTFR
write(1,*) "flux("//trim(chartu)//")^1), error, relative error (%)"
write(1,'(i1,3f16.6)') 1, flux1,flux1*relerrflux1,100*relerrflux1
write(1,'(i1,3f16.6)') 2, flux2,flux2*relerrflux2,100*relerrflux2
write(1,'(i1,3f16.6)') 3, flux3,flux3*relerrflux3,100*relerrflux3
write(1,*) "MDcycles:",MDcycles
write(1,*) "NSA,NSB,NOSA,NOSB",NSA,NSB,NOSA,NOSB
write(1,*) "NCROSS",NCROSS(1:3)
write(1,*) "NEFFCROSS",NEFFCROSS(1:3)
write(1,*) "NEFFCROSS/NCROSS",1.d0*NEFFCROSS(1:3)/NCROSS(1:3)
if ((NOSA /=MDcycles).AND.(NSB/=MDcycles)) then
  write(1,*) "reactive flux correction: c1=NOSA/(MDcycles-NSB)"
  c1=1.d0*NOSA/(MDcycles-NSB)
  write(1,*) "c1=",c1
  write(1,*) "flux1*c1,flux2*c1,flux3*c1:"
  write(1,'(3f16.6)')flux1*c1,flux2*c1,flux3*c1
  if (POTENTIAL=="HARMOSC") then
    write(1,*) "units kappa might need a check!!!"
    kapTST=1/sqrt(2*pi*mass_eVns2_A2*beta)
    write(1,*)  "kappa TST result", kapTST 
    PxinAL=exp(-beta*.5*kharm*MD_INTFL**2)/&
    (sqrt( pi/(2*beta*kharm) ) * erfc01( - sqrt( (beta*kharm)/2 )*MD_INTFR  ))
    PxinAM=exp(-beta*.5*kharm*MD_INTFM**2)/&
    (sqrt( pi/(2*beta*kharm) ) * erfc01( - sqrt( (beta*kharm)/2 )*MD_INTFR  ))
    PxinAR=exp(-beta*.5*kharm*MD_INTFR**2)/&
    (sqrt( pi/(2*beta*kharm) ) * erfc01( - sqrt( (beta*kharm)/2 )*MD_INTFR  ))
    write(1,'(a22,3f16.6)') "PxinAL,PxinAM,PxinAR",PxinAL,PxinAM,PxinAR
    write(1,*) "TST results:"
    write(1,'(3f16.6)') kapTST*PxinAL,kapTST*PxinAM,kapTST*PxinAR 
  endif
endif
pMD1=flux1*dt_ut
pMD2=flux2*dt_ut
pMD3=flux3*dt_ut
write(1,*) "flux1: one effective crossing per", 1.d0/pMD1, " timesteps"
write(1,*) "flux2: one effective crossing per", 1.d0/pMD2, " timesteps"
write(1,*) "flux3: one effective crossing per", 1.d0/pMD3, " timesteps"
write(1,*) "pMD1",pMD1, "(1-p)/p",(1.d0-pMD1)/pMD1
write(1,*) "pMD2",pMD2, "(1-p)/p",(1.d0-pMD2)/pMD2
write(1,*) "pMD3",pMD3, "(1-p)/p",(1.d0-pMD3)/pMD3
teffMD1=relerrflux1**2*MDcycles
teffMD2=relerrflux2**2*MDcycles
teffMD3=relerrflux3**2*MDcycles
write(1,*) "efficiency time 1",teffMD1
write(1,*) "efficiency time 2",teffMD2
write(1,*) "efficiency time 3",teffMD3
corMd1=teffMD1*pMD1/(1.d0-pMD1)
corMd2=teffMD2*pMD2/(1.d0-pMD2)
corMd3=teffMD3*pMD3/(1.d0-pMD3)
write(1,*) "correlation 1",corMd1
write(1,*) "correlation 2",corMd2
write(1,*) "correlation 3",corMd3
write(1,*) "\end{verbatim}"
write(1,*) "\clearpage"
write(1,*) "\begin{verbatim}"

write(1,*) "transmission results: !! units might need a check!!!!!!!!!"
write(1,*) "trans_INTFL=",trans_INTFL
write(1,*) "trans_INTFM=",trans_INTFM
write(1,*) "trans_INTFR=",trans_INTFR
write(1,*) "(unnormalized) kappa (A/ns), error, relative error (%)"
write(1,'(3f16.6)') kappa,rekap*kappa,rekap*100 
print *,"~~~~~~~~~~~~~~~"
print *,"kappa (A/ns), error, %error", kappa,rekap*kappa,rekap*100
write(1,*) "(unnoormalized) kappaTST (A/ns), error, relative error (%)"
write(1,'(3f16.6)') kappaTST,rekapTST*kappaTST,rekapTST*100
kapTSTtheor=1/sqrt(2*pi*mass_eVns2_A2*beta)
write(1,*) "exact TST result", kapTSTtheor 

!??kappaTSTBrownian=sqrt(1.d0/(pi*mass_eVns2_A2*beta*gamma*dt))
!??? what is this formula?? 
! here the formula from Bernd's thesis
write(1,*) "normalized kappa (k/kTST):", kappa/kapTSTtheor
print *,"normalized kappa", kappa/kapTSTtheor
print *,"PATS:",prob1

if (POTENTIAL=="DOUBLEWELL") then
  wb=-DOUBLEWELLk2
  kappaTSTBrownian=(1.d0/wb)*(-gamma/2.d0+sqrt(gamma**2/4.d0+wb**2))
  write(1,*) "exact Brownian dynamics TST result", kappaTSTBrownian
endif

write(1,*) "num. integ. P_A(TS)=",prob1,prob2
if (POTENTIAL=="HARMOSC") then
  PxinATS=exp(-beta*.5*kharm*trans_INTFM**2)/&
  (sqrt( pi/(2*beta*kharm) ) * erfc01( - sqrt( (beta*kharm)/2 )*trans_INTFM  ))
  PxinATS2=exp(-beta*.5*kharm*trans_INTFM**2)/&
  (sqrt( pi/(2*beta*kharm) ) * erfc01( - sqrt( (beta*kharm)/2 )*trans_INTFR  ))
  write(1,*) "P_A(TS)=",PxinATS,PxinATS2
else
  PxinATS=prob1
  PxinATS2=prob2
endif

if (POTENTIAL=="DOUBLEWELL") then
  write(1,*) "results numerical integration double-well"
  tol=1.d-9
  PxinATS=0.5d0 !on to only halve a bin
  x=0.d0
  do
    x=x+dxinteg
    prob1=exp(-beta*(doublewellk4*x**4+doublewellk2*x**2))
    if (prob1<tol) exit
    PxinATS=PxinATS+prob1
  enddo
  PxinATS=PxinATS*dxinteg
  PxinATS=1.d0/PxinATS
  write(1,*)"probability to be at x=0 given x<=0 equals:",PxinATS 
endif
 
write(1,*) "reactive flux result rate"
write(1,*) kappa*PxinATS,kappa*PxinATS2
print *,"k=",kappa*PxinATS
if (two_point_method) then
  write(1,*) "two-point method:"
  if (brownian) then
    write(1,*) kappa*PxinATS*kappaTSTBrownian,kappa*PxinATS2*kappaTSTBrownian
  else
    write(1,*) kappa*PxinATS*kapTST,kappa*PxinATS2*kapTST
  endif
endif
print *,"~~~~~~~~~~~"
write(1,*) "\end{verbatim}"
write(1,*) "\clearpage"
write(1,*) "\begin{verbatim}"



if (numint==0) return
write(1,*) "TIS results"
totrelerrPcross=0
totsimPcross=0
totPcross=1
teffPcrossOpt=0.d0
do i=1,numint
  write(1,*) "ensemble",i
  idetect=intfr(i);if (i<numint) idetect=intfm(i+1)
  write(1,'(a16,4f10.4)') "IL,IM,IR,Idetect",intfl(i),intfm(i),intfr(i),idetect
  write(1,'(a20,3f16.6)') "pcross,err,rel. e(%)",crossprobab(i), &
  crossprobab(i)*relerrcr(i), 100*relerrcr(i)
  write(1,'(a20,3f16.6)') "Lacc,Ltr, Ltr/Lacc",Lacc(i),Ltr(i), Ltr(i)/Lacc(i)
  TISsim(i)=TIScycles(i)*Ltr(i)
  teffTIS(i)=relerrcr(i)**2*TISsim(i)
  TisCyclesOpt(i)=sqrt(teffTIS(i))/Ltr(i)
  write(1,*) "TIScycles, Totsim",TIScycles(i),TISsim(i)
  write(1,*) "accRate, cor",accrate(i),corcr(i)
  write(1,*) "teff", teffTIS(i)
  totpcross=totpcross*crossprobab(i)
  totrelerrPcross= totrelerrPcross+relerrcr(i)**2
  totsimPcross=totsimPcross+TISsim(i)
  teffPcrossOpt=teffPcrossOpt+sqrt(teffTIS(i))
enddo
totrelerrPcross=sqrt(totrelerrPcross)
teffPcrossOpt=teffPcrossOpt**2

if (TASK=="PPS") then
  write(1,*) "\end{verbatim}"
  write(1,*) "\clearpage"
  write(1,*) "\begin{verbatim}"
  write(1,*) "additional PPS 000/001 results:"
  write(1,'(a20,3f16.6)') "T0,err,rel.e(%)",T0,T0*reT0,100*reT0
  write(1,'(a20,3f16.6)') "T1,err,rel.e(%)",T1,T1*reT1,100*reT1
  write(1,'(a20,3f16.6)') "LtrT0, Ltr/Lacc",LTRT0,LtrT0/(T0+2)
  write(1,'(a20,3f16.6)') "LtrT1, Ltr/Lacc",LTRT1,LtrT1/(T1+2)
  simT0=NL0*LtrT0
  simT1=NL1*LtrT1
  teffT0=reT0**2*simT0
  teffT1=reT1**2*simT1
  CyclesOptT0=sqrt(teffT0)/LtrT0
  write(1,*) "cyclesT0, Totsim",NL0,simT0
  write(1,*) "cyclesT1, Totsim",NL1,simT1
  write(1,*) "accRateT0, cor",accrateT0,NcorT0
  write(1,*) "accRateT1, cor",accrateT1,NcorT1
  write(1,*) "teffT0", teffT0
  write(1,*) "teffT1", teffT1
  fluxPPS=1.d0/((T0+T1)*dt_ut)
  rePPS=(reT0*T0+reT1*T1)/(T0+T1)
  teffPPS=rePPS**2*simT0
endif


write(1,*) "\end{verbatim}"
write(1,*) "\clearpage"
write(1,*) "\begin{verbatim}"
write(1,*) "Combining all results:"
if ( (task/="PPS").AND.(numint>0).AND.(MD_INTFL==intfl(1)).AND.&
                                      (MD_INTFL==intfm(1)) ) then
  flux=flux1
  re_flux=relerrflux1
  teffflux=teffMD1
  fluxcycles=MDcycles
  totsimflux=1.d0*fluxcycles
  CyclesOptflux=sqrt(teffflux)
else if ((task/="PPS").AND. (numint>0).AND.(MD_INTFL==intfl(1)).AND.&
                                           (MD_INTFM==intfm(1)) ) then
  flux=flux2
  re_flux=relerrflux2
  teffflux=teffMD2
  fluxcycles=MDcycles
  totsimflux=1.d0*fluxcycles
  CyclesOptflux=sqrt(teffflux)

else if (TASK=="PPS") then
  flux=fluxPPS
  re_flux=rePPS
  teffflux=teffPPS
  fluxcycles=NL0
  totsimflux=simT0
  CyclesOptflux=sqrt(teffflux)/LTRT0
else
 print *,"flux doesnot match TIS simulations" 
endif

write(1,*) "flux=",flux,"+-",100*re_flux,"%"
write(1,*) "Pcross=",totpcross,"+-",100*totrelerrPcross,"%"
rate=flux*totpcross
relerrRate=sqrt(re_flux**2+totrelerrPcross**2)
write(1,*) "rate=",rate,"+-",100*relerrRate,"%"
write(1,*) "---------------------------------------"
totsimRate=totsimflux+totsimPcross
teffPcross=totrelerrPcross**2*totsimPcross
teffRate=relerrRate**2*totsimRate
write(1,*) "flux:sim.time, teff",totsimflux,teffflux
write(1,*) "Pcross:sim.time, teff",totsimPcross,teffPcross
write(1,*) "Rate:sim.time, teff",totsimRate,teffRate
write(1,*) "Optimeze Pcross teff:",teffPcrossOpt
teffRateOpt=(sqrt(teffflux)+sqrt(teffPcrossOpt))**2
write(1,*) "Optimeze rate teff:",teffrateOpt
write(1,*) "ratio of chosen NCYCLES"
Ncycmin=minval(TIScycles)
Ncycmin=min(1.d0*fluxcycles,Ncycmin)
write(1,*) "flux:",1.d0*fluxcycles/Ncycmin
write(1,*) "TIS:"
do i=1,numint
  write(1,*) i,1.d0*TIScycles(i)/Ncycmin
enddo
cyclesOPTmin=minval(TIScyclesOpt)
cyclesOPTmin=min(cyclesOPTmin,cyclesoptflux)
write(1,*) "ratio of optimal NCYCLES"
write(1,*) "flux:",cyclesoptflux/cyclesOPTmin
write(1,*) "TIS:"
do i=1,numint
  write(1,*) i,TIScyclesOpt(i)/cyclesOPTmin
enddo


end subroutine write_results
!ES-------------------------------------------------------------

!EM----------------------------------------------------------------------
end Module summary

