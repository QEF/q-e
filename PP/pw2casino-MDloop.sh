#!/bin/bash

###################################################################################################
#
# (C) 2010 by Norbert Nemec <Norbert@Nemec-online.de>
#
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
###################################################################################################
#
# Patch script to write a series of wave functions in CASINO format
# at the consecutive time steps during a MD run.
#
# To implement this functionality, the output routine needs to be called from the main MD loop.
# This would erode the current separation of functionality in PW/ and in PP/
# As the functionality is also rather specialized, this script is provided to apply the
# necessary changes in an existing distribution.
#
# USAGE:
#
# * Call the script a
#     ./pw2casino-MDloop.sh
#   from within the PP/ directory.
#
# * Check for potential error messages. Errors are most likely due to changes in the distribution.
#   In that case, this script needs to be updated. Please contact the author for help if necessary.
#
# * Compile the program as usual.
#
# Beware that this patch does not provide any run-time options for simplicity. Edit the arguments
# in the function call to write_casino_wfn in the PW/pwscf.f90 file if necessary.
#
###################################################################################################

PPDIR=$(dirname $0)
SCRIPTNAME=$(basename $0)
cd $PPDIR/..

PATCHFILE=pw2casino-MDloop.diff

cat > $PATCHFILE <<'EOF'
patch files to call pw2casino from MD loop

From: Mike Towler <mdt26@cam.ac.uk>


---
 PP/Makefile  |    4 ++--
 PW/Makefile  |    2 ++
 PW/pwscf.f90 |   15 +++++++++++++++
 3 files changed, 19 insertions(+), 2 deletions(-)

diff --git a/PP/Makefile b/PP/Makefile
index 895b70e..0456149 100644
--- a/PP/Makefile
+++ b/PP/Makefile
@@ -120,9 +120,9 @@ projwfc.x : projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS)
 		projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
 	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
-pw2casino.x : pw2casino.o pw2casino_write.o pw2blip.o $(PPOBJS) $(MODULES) $(LIBOBJS)
+pw2casino.x : pw2casino.o $(PPOBJS) $(MODULES) $(LIBOBJS)
 	$(LD) $(LDFLAGS) -o $@ \
-		pw2casino.o pw2casino_write.o pw2blip.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
+		pw2casino.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
 	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
 pw2wannier90.x : pw2wannier90.o $(PPOBJS) $(MODULES) $(LIBOBJS)
diff --git a/PW/Makefile b/PW/Makefile
index bc7d671..4a4ee64 100644
--- a/PW/Makefile
+++ b/PW/Makefile
@@ -141,6 +141,7 @@ print_ks_energies.o \
 punch.o \
 pw_restart.o \
 pwcom.o \
+pw2blip.o \
 qvan2.o \
 rcgdiagg.o \
 rdiagh.o \
@@ -224,6 +225,7 @@ wannier_check.o \
 wannier_clean.o \
 wannier_occ.o \
 wannier_enrg.o \
+write_casino_wfn.o 
 
 QEMODS=../Modules/libqemod.a
 
diff --git a/PW/pwscf.f90 b/PW/pwscf.f90
index bb923bd..7418e21 100644
--- a/PW/pwscf.f90
+++ b/PW/pwscf.f90
@@ -29,6 +29,8 @@ PROGRAM pwscf
   USE open_close_input_file_interf, ONLY : open_input_file, close_input_file
   !
   IMPLICIT NONE
+  CHARACTER(4) postfix
+  INTEGER itercount
   !
   LOGICAL :: xmlinput = .false.
   CHARACTER (len=iotk_attlenx) :: attr
@@ -91,8 +93,11 @@ PROGRAM pwscf
   !
   CALL init_run()
   !
+  itercount=0
+  !
   main_loop: DO
      !
+     itercount=itercount+1
      ! ... electronic self-consistentcy
      !
      CALL electrons()
@@ -102,6 +107,16 @@ PROGRAM pwscf
        CALL stop_run( conv_elec )
      ENDIF
      !
+     write(postfix,'(i4.4)')itercount
+     CALL write_casino_wfn( &
+     .false., & ! gather
+     .true.,  & ! blip
+     1.0d0,   & ! multiplicity
+     .true.,  & ! binwrite
+     .true.,  & ! single_precision_blips
+     0,       & ! n_points_for_test
+     '.'//postfix)   ! postfix
+     !
      ! ... if requested ions are moved
      !
 !  CALL ions()
EOF

patch -p 1 --dry-run -f -p 1 < $PATCHFILE

if [ "$?" != "0" ] ; then
    echo
    echo "ERROR: Patch did not apply cleanly."
    echo "--> Very likely, the distribution has been changed without updating the script PP/$SCRIPTNAME"
    echo "--> Either apply the patch $PATCHFILE manually (it is fairly simple) or contact the author for help:"
    echo "-->     <mdt26@cam.ac.uk>"
    exit 1
fi

patch -p 1 -s -b < $PATCHFILE
rm -f $PATCHFILE

mv PP/pw2blip.f90 PW/pw2blip.f90
mv PP/pw2casino_write.f90 PW/write_casino_wfn.f90
install/makedeps.sh
