#Debugging flags for the Intel compilers:
#-g 	Produces symbolic debug information in the object files.
#-CB 	Performs array bounds checking.
#-u 	For Fortran only: same as setting IMPLICIT NONE throughout the code.
#-traceback 	For Fortran only: enables display of source file traceback information at run time when a severe error occurs.
#Optimization flags
#-O0  	 Disables optimization
#-O1 	Enables optimizations for speed, while being aware of code size.
#-O2 or -O 	Enables optimizations.
#-O3 	Enables -O2 plus more aggressive optimizations. Use with care.
#-i8 	Sets default size of INTEGER to 8 bytes.
#-r8 	Sets default size of REAL to 8 bytes.


#F90 = f90 
#F90= ifort
F90=gfortran
#FFLAGS=    -g -CB -u -traceback
.SUFFIXES: .f90 .f

OBJS1 = stringlengths.o defaultfile.o types.o assign_objects.o alloc_objects.o inputpar.o object_modules.o random.o charext.o shell.o read_input_file.o arimetric.o inquire.o restart.o cpmd_subr.o parameters.o atoms.o initharm.o initpbd.o initwca.o initexternal.o  initnosehoover.o initlange.o initdynext.o  initnuminteg.o
OBJS2 = forcefieldharm.o forcefieldnicolis.o forcefieldpbd.o forcefieldwca.o forcefieldiot.o forcefield.o  orderparameter.o thermo.o mdintegrators.o mdstep.o randomnize.o convert.o initmd.o  initiot.o  initpot.o 
OBJS3 = wannier.o outmd.o  outtis.o outpps.o propagator.o biasvacc.o tisstep.o transstep.o null.o swapfm.o swap.o ppsstep.o
OBJS4 = inittis.o initpps.o init.o  interactive.o md.o tis.o pps.o trans.o ffs.o numericinteg.o exec.o close.o main.o
OBJS = $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
PRGM = TISMOL.exe 

$(PRGM):  $(OBJS)
	$(F90) $(FFLAGS) -o $(PRGM) $(OBJS)


.f90.o:
	$(F90) $(FFLAGS) -c $<
.f.o  :
	$(F77) $(FFLAGS) -c $<

clean:
	rm *.o *.mod *.dat *.xyz fort.* 
cleanall:
	rm *.o *.mod *.dat *.xyz fort.* 0*/*dat 0*/*xyz 0*/input.TISMOL md/*dat md/*xyz *exe md/input.TISMOL
