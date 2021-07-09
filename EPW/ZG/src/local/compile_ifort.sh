ifort -g -c JDOS_Gaus.f90
ifort -o JDOS_Gaus.x JDOS_Gaus.o
#
ifort -g -c create_qlist.f90
ifort -o create_qlist.x create_qlist.o
#
ifort -g -c rotate.f90
ifort -o rotate.x rotate.o
#
ifort -g -c kpoints_band_str_unfold.f90
ifort -o kpoints_band_str_unfold.x kpoints_band_str_unfold.o
