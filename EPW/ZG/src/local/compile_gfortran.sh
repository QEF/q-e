gfortran -g -c JDOS_Gaus.f90
gfortran -o JDOS_Gaus.x JDOS_Gaus.o
#
gfortran -g -c create_qlist.f90
gfortran -o create_qlist.x create_qlist.o
#
gfortran -g -c rotate.f90
gfortran -o rotate.x rotate.o
#
gfortran -g -c kpoints_band_str_unfold.f90
gfortran -o kpoints_band_str_unfold.x kpoints_band_str_unfold.o
