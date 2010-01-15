# Makefile for PW
#Modified for v4.0 by obm


include ../../make.sys
IFLAGS         = -I../../include
LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../Multigrid/mglib.a



# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../EE $(MOD_FLAG) . $(MOD_FLAG)../../PW \
          $(MOD_FLAG)../../PH

PWOBJS = ../../PW/libpw.a
                              
EEOBJS = \
../../EE/dvdr_tao.o                           \
../../EE/add_dcc_field.o                      \
../../EE/v_h_from_rho_r.o                     \
../../EE/multiscale.o                         \
../../EE/writetofile.o                        \
../../EE/mg_pb_solver.o                       \
../../EE/add_boundary.o                       \
../../EE/add_dccdil_forces.o                  \
../../EE/init_ee.o                            \
../../EE/calc_ecomp.o                         \
../../EE/setlocalcoul.o                       \
../../EE/set_fft_dim_coarse.o                 \
../../EE/set_mltgrid_dim.o                    \
../../EE/data_structure_coarse.o              \
../../EE/add_ele_corr.o                       \
../../EE/write_ee_summary.o                   \
../../EE/ggen_coarse.o

PHOBJS=../../PH/libph.a

                                                                
                                     
EEMODS=../../EE/ee_mod.o ../../EE/gcoarse_mod.o

LROBJS= \
lr_variables.o \
lr_charg_resp.o \
lr_newd.o \
bcast_lr_input.o \
lr_readin.o \
lr_alloc_init.o \
lr_calc_dens.o \
lr_dot.o \
lr_dealloc.o \
lr_ortho.o \
lr_read_wf.o \
lr_normalise.o \
lr_lanczos.o \
lr_apply_liouvillian.o \
lr_main.o \
lr_dv_setup.o \
lr_setup_dgc.o \
lr_solve_e.o \
lr_dvpsi_e.o \
lr_ch_psi_all.o \
lr_cgsolve_all.o \
lr_h_psiq.o \
lr_sm1_psi.o \
stop_lr.o \
lr_read_d0psi.o \
lr_restart.o \
lr_write_restart.o \
print_clock_lr.o \
sd0psi.o \
lr_set_boxes_density.o \
lr_init_nfo.o 
 
MODULES = ../../Modules/*.o

all :  tddfpt.x

#tddfpt.x : $(PWOBJS) $(EEOBJS) $(PHOBJS) $(LROBJS) 
#	$(MPIF90) $(LDFLAGS) -o tddfpt.x \
#	$(PWOBJS) $(EEOBJS) $(PHOBJS) $(LROBJS) $(MODULES) $(EEMODS) $(LIBOBJS) $(LIBS) 
#	- ( cd ../../bin; ln -fs ../../TDDFPT/tddfpt.x . )
#
tddfpt.x : $(LROBJS) $(PWOBJS) $(MODULES) $(PHOBJS) $(PWOBJS) $(EEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		$(LROBJS)  $(MODULES) $(PWOBJS) $(PHOBJS) $(EEMODS) $(LIBOBJS) $(LIBS)

dene.x : obm_tests.o 
	$(MPIF90) $(LDFLAGS) -o dene.x \
        $(PWOBJS) $(PHOBJS) obm_tests.o $(MODULES) $(LIBOBJS) $(LIBS)
clean :
	- /bin/rm -f tddfpt.x *.o *~ *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
