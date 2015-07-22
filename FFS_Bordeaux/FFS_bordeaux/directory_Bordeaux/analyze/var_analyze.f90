!BM---------------------------------------------
Module var_analyze
implicit none
character*100::finpa,fileMD,fordp,filecross,fcross2,fpath,fpacc,fexc,fsum
character*100, allocatable::PSFILES(:)
logical::MDFILE,pathfile,ORDFILE,crossFILE
integer::iff,NFMAX,NLinpa,NLMD,NLord,NLpath
integer::NLcross,ngrid,Lpmax,lbmax,Nc2,NFpath,Nlp
double precision::pi, unit_t,ETav
double precision, allocatable::cross(:),actEN(:),ordpath(:)
double precision::FLUX,Pctot,rate,Ncfl,stdfl,rhoHA,rhoHB
integer, allocatable::apcro1(:),NLFpath(:),skp(:),NLPA(:)
double precision, allocatable:: lamA(:),lamB(:),lamC(:)
integer:: NHA,Ncross
double precision::errHAM, errTemp, errflux,rerrflux,rerrk,errk,rerrpc,errpc
character*3, allocatable::EXT(:)
character*100, allocatable::ficr(:),fiactEN(:)
double precision, allocatable::pcross(:),errpcr(:),rerrpcr(:)
double precision, allocatable::stdpc(:),Ncpc(:),LMC(:),maxp(:),minp(:)
double precision, allocatable::arrl(:),arrPc(:),Neff(:)
character*100::dir
integer::NLpctot
end Module var_analyze
!EM---------------------------------------------

