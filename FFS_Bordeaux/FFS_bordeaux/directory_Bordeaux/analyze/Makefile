#F90 = f90 
#F90 = ifort
F90 = gfortran
FFLAGS=   
.SUFFIXES: .f90 .f


OBJS1 = var.o erfc.o shell.o read_input_file.o param.o inquire.o init.o gnuplot.o data_analysis.o funcplot.o output.o GAMMA_func.o analyze_EN.o analyze_OP.o analyze_pcross.o analyze_flux.o MDanalyze.o transmission.o numintanalyze.o analyze_shoot.o Pathanalyze.o match.o summary.o movies.o main.o

OBJS = $(OBJS1) 
PRGM = analyze.exe 

$(PRGM):  $(OBJS)
	$(F90) $(FFLAGS) -o $(PRGM) $(OBJS)


.f90.o:
	$(F90) $(FFLAGS) -c $<
.f.o  :
	$(F77) $(FFLAGS) -c $<

clean:
	rm *.o *.mod *dat *tmp *gnu *eps fort.*  
