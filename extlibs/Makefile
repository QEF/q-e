
#===============================
# espresso
#===============================
#
include ../make.sys

#
# MAIN target
#
all: 
	$(MAKE) libblas_$(BLAS_LIBS_SWITCH) 
	$(MAKE) liblapack_$(LAPACK_LIBS_SWITCH) 


# BLAS
libblas : libblas_$(BLAS_LIBS_SWITCH)

libblas_external : fake_blas_external

fake_blas_external :
	touch fake_libblas.a
	rm fake_libblas.a

libblas_internal:
	if test ! -d ../BLAS; then \
	( tar -xvf archive/blas-1.tar -C ../); fi
	if test -e install/make_blas.inc; then \
	(cp install/make_blas.inc ../BLAS/make.inc; \
	cd ../BLAS; $(MAKE) $(BLASLIB)); else \
        (echo "no configuration file found for blas"; \
	echo "run configure from main QE dir"); fi

# LAPACK
liblapack : liblapack_$(LAPACK_LIBS_SWITCH) 

liblapack_external : fake_lapack_external

fake_lapack_external :
	touch fake_liblapack.a
	rm fake_liblapack.a

liblapack_internal:
	if test ! -d ../lapack-3.2; then \
	( tar -xvf archive/lapack-3.2.tar -C ../ ); fi
	if test -e install/make_lapack.inc; then \
	(cp install/make_lapack.inc ../lapack-3.2/make.inc; \
	cd ../lapack-3.2; $(MAKE) lapacklib); else \
        (echo "no configuration file found for lapack"; \
        echo "run configure from main QE dir"); fi

#Multigrid
mglib:
	if (test ! -d ../Multigrid || \
	test ! -e ../Multigrid/Makefile); then \
	( tar -xvf archive/multigrid.tar -C ../ ; \
	cd ../Multigrid; $(MAKE) all); fi
###################################
# cleaning
###################################
# each lib independently
blas_clean:
	if test -d ../BLAS; then (cd ../BLAS; $(MAKE) clean); fi
blas_veryclean:
	if test -d ../BLAS; then (rm -R -f ../BLAS); fi
lapack_clean:
	if test -d ../lapack-3.2; then (cd ../lapack-3.2; $(MAKE) clean); fi
lapack_veryclean:
	if test -d ../lapack-3.2; then (rm -R -f ../lapack-3.2); fi
multigrid_clean:
	if test -d ../Multigrid; then (cd ../Multigrid; $(MAKE) clean); fi
multigrid_veryclean:
	if test -d ../Multigrid; then (rm -R -f ../Multigrid); fi

# general cleaning
clean:
	if test -d ../lapack-3.2 ;  then ( cd ../lapack-3.2;    $(MAKE) $@ ) ; fi
	if test -d ../BLAS ;  then ( cd ../BLAS;    $(MAKE) $@ ) ; fi
	if test -d ../Multigrid ;  then ( cd ../Multigrid;    $(MAKE) $@ ) ; fi

veryclean:
	if test -d ../BLAS; then ( rm -R -f ../BLAS); fi
	if test -d ../lapack-3.2 ; then ( rm -R -f ../lapack-3.2); fi
	if test -d ../Multigrid ; then ( rm -R -f ../Multigrid); fi
