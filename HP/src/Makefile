# Makefile for HP

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
	  $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules
HPOBJS = \
hp_allocate_q.o \
hp_bcast_input.o \
hpcom.o \
hp_clean_q.o \
hp_close_q.o \
hp_calc_chi.o \
hp_check_pert.o \
hp_check_type.o \
hp_dnstot_sum_q.o \
hp_dnsq.o \
hp_dealloc_q.o \
hp_dealloc_1.o \
hp_dealloc_2.o \
hp_dvpsi_pert.o \
hp_efermi_shift.o \
hp_find_inequiv_sites.o \
hp_generate_grids.o \
hp_init.o \
hp_init_q.o \
hp_load_q.o \
hp_ns_trace.o \
hp_openfil_q.o \
hp_psymdvscf.o \
hp_postproc.o \
hp_prepare_q.o \
hp_print_clock.o \
hp_q_points.o \
hp_R_points.o \
hp_rotate_dnsq.o \
hp_readin.o \
hp_read_chi.o \
hp_read_dnsq.o \
hp_run_nscf.o \
hp_solve_linear_system.o \
hp_symdnsq.o \
hp_symdvscf.o \
hp_setup_q.o \
hp_summary.o \
hp_summary_q.o \
hp_stop_smoothly.o \
hp_write_chi.o \
hp_write_chi_full.o \
hp_write_dnsq.o 

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

TLDEPS= hplibs

all : tldeps libs-hp hp.x 

libs-hp : libhp.a 

hp.x : hp_main.o libhp.a $(PWOBJS) $(LRMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ hp_main.o libhp.a \
		$(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../HP/src/hp.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libhp.a : $(HPOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/hp.x 

include make.depend
# DO NOT DELETE
