!BM------------------------------------
module atoms
implicit none

CONTAINS

!BS----------------------------------------
subroutine set_atomic_masses(masses,atom_types,Npart)
implicit none
integer, intent(in)::Npart
character*2, intent(in)::atom_types(Npart)
double precision, intent(out)::masses(Npart)
character*2::Periodic_system(100)
double precision::atommasses(100)
integer::i,j

print *,"DEFINING MASS VALUES FROM PERIODIC SYSTEM"

atommasses=(/ 1.00797D0,  4.0026D0,    6.939D0,   9.0122D0,  10.811D0, &
             12.01115D0, 14.0067D0,   15.9994D0, 18.9984D0,  20.183D0, &
             22.9898D0,  24.312D0,    26.9815D0, 28.086D0,   30.9738D0,&
             32.064D0,   35.453D0,   39.948D0,   39.102D0,   40.080D0, &
             44.956D0,   47.900D0,   50.942D0,   51.996D0,   54.938D0, &
             55.847D0,   58.933D0,   58.710D0,   63.540D0,   65.370D0, &
             69.720D0,   72.590D0,   74.922D0,   78.960D0,   79.909D0, &
             83.800D0,   85.470D0,   87.620D0,   88.905D0,   91.220D0, &
             92.906D0,   95.940D0,   98.000D0,  101.070D0,  102.905D0, &
            106.400D0,  107.870D0,  112.400D0,  114.820D0,  118.690D0, &
            121.750D0,  127.600D0,  126.904D0,  131.300D0,  132.905D0, &
            137.340D0,  138.910D0,  140.120D0,  140.907D0,  144.240D0, &
            147.000D0,  150.350D0,  151.960D0,  157.250D0,  158.924D0, &
            162.500D0,  164.930D0,  167.260D0,  168.934D0,  173.040D0, &
            174.970D0,  178.490D0,  180.948D0,  183.850D0,  186.200D0, &
            190.200D0,  192.200D0,  195.090D0,  196.967D0,  200.590D0, &
            204.370D0,  207.190D0,  208.980D0,  210.000D0,  210.000D0, &
            222.000D0,  250.000D0,  250.000D0,  250.000D0,  250.000D0, &
            250.000D0,  250.000D0,  250.000D0,  250.000D0,  250.000D0, &
            250.000D0,  250.000D0,  250.000D0,  250.000D0,  250.000D0 /)

Periodic_system=(/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne",    &
                  "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca",    &
                  "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn",    &
                  "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr",    &
                  "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",    &
                  "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd",    &
                  "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",    &
                  "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg",    &
                  "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",    &
                  "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm"    /)    

do i=1,Npart
  do j=1,100
    if (trim(atom_types(i))==trim(Periodic_system(j))) then
      masses(i)=atommasses(j)
      print *,i,atom_types(i),masses(i)
      exit 
    endif
    if (j==100) then
      print *,"ERROR set_atomic_masses atom_types(i)=",atom_types(i)
      stop
    endif
  enddo
enddo

end subroutine set_atomic_masses
!ES----------------------------------------

end module atoms
!EM-------------------------------------
