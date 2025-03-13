.SUFFIXES : .o .f90 .F90

%.o: %.f90
	$(FXX) $(FXXOPT) -c $<
%.o: %.F90
	$(FXX) $(FXXOPT)  -c $<

LIBMBD_C_API = 0
OBJS := mbd.o mbd_c_api.o mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_defaults.o mbd_density.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_methods.o mbd_rpa.o mbd_scs.o mbd_ts.o mbd_utils.o mbd_version.o mbd_vdw_param.o
ifeq ($(LIBMBD_C_API),0)
OBJS := $(filter-out mbd_c_api.o,$(OBJS))
endif

libmbd.a: $(OBJS)
	ar -r $@ $^


mbd.o: mbd_constants.o mbd_damping.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_methods.o mbd_ts.o mbd_utils.o mbd_vdw_param.o mbd_defaults.o mbd_version.o
mbd_c_api.o: mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_methods.o mbd_ts.o mbd_utils.o mbd_density.o mbd_version.o
mbd_coulomb.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_lapack.o mbd_linalg.o mbd_matrix.o
mbd_damping.o: mbd_constants.o mbd_gradients.o mbd_utils.o mbd_defaults.o
mbd_defaults.o: mbd_constants.o
mbd_density.o: mbd_constants.o mbd_geom.o mbd_linalg.o mbd_lapack.o
mbd_dipole.o: mbd_constants.o mbd_damping.o mbd_geom.o mbd_gradients.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_utils.o
mbd_formulas.o: mbd_constants.o mbd_gradients.o mbd_utils.o
mbd_geom.o: mbd_constants.o mbd_lapack.o mbd_utils.o mbd_vdw_param.o mbd_defaults.o
mbd_gradients.o: mbd_constants.o
mbd_hamiltonian.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_utils.o
mbd_lapack.o: mbd_constants.o mbd_utils.o
mbd_linalg.o: mbd_constants.o
mbd_matrix.o: mbd_constants.o mbd_lapack.o mbd_utils.o
mbd_methods.o: mbd_constants.o mbd_damping.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_rpa.o mbd_scs.o mbd_utils.o
mbd_rpa.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_matrix.o mbd_utils.o
mbd_scs.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_utils.o
mbd_ts.o: mbd_constants.o mbd_damping.o mbd_geom.o mbd_utils.o
mbd_utils.o: mbd_constants.o mbd_gradients.o
mbd_vdw_param.o: mbd_constants.o mbd_utils.o

mbd_version.o: mbd_version.f90 mbd_constants.o
	$(FXX) $(FXXOPT) -c $<

mbd_version.f90: mbd_version.f90.in
	sed  -e 's/@PROJECT_VERSION_MAJOR@/0/' mbd_version.f90.in  > mbd_version.f90
	sed -i -e 's/@PROJECT_VERSION_MINOR@/13/'  mbd_version.f90
	sed -i -e 's/@PROJECT_VERSION_PATCH@/0/' mbd_version.f90
	sed -i -e 's/@PROJECT_VERSION_SUFFIX@/0/' mbd_version.f90

.PHONY: clean distclean
clean:
	rm -f *.a *.o mbd_version.f90

distclean: clean
	rm -f *.mod
	rm -f $(LIB)
