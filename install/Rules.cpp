.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	$(CPP) $(CPPFLAGS) $*.f90 $*.F90
	$(F90) $(F90FLAGS) -c $*.F90 -o $*.o

.f.o:
	$(F77) $(F77FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<
