.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	$(CPP) $(CPPFLAGS) $*.f90 $*.F90
	$(F90) $(F90FLAGS) $(MODULEFLAG) -c $*.F90 -o $*.o
#	rm $*.F90

.f.o:
	$(F77) $(F77FLAGS) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<
