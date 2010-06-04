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

From: Norbert Nemec <Norbert@Nemec-online.de>


---
 PP/Makefile  |    4 ++--
 PW/Makefile  |    4 +++-
 PW/pwscf.f90 |   15 +++++++++++++++
 3 files changed, 20 insertions(+), 3 deletions(-)

diff --git a/PP/Makefile b/PP/Makefile
index ec24bc6..29dbceb 100644
--- a/PP/Makefile
+++ b/PP/Makefile
@@ -123,9 +123,9 @@ projwfc.x : projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS)
 		projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
 	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
-pw2casino.x : pw2casino.o pw2casino_write.o pw2blip.o $(PPOBJS) $(MODULES) $(LIBOBJS)
+pw2casino.x : pw2casino.o $(PPOBJS) $(PWOBJS) $(MODULES) $(EEMODS) $(LIBOBJS)
 	$(LD) $(LDFLAGS) -o $@ \
-		pw2casino.o pw2casino_write.o pw2blip.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
+		pw2casino.o $(PPOBJS) $(MODULES) $(EEMODS) $(PWOBJS) $(LIBOBJS) $(LIBS)
 	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
 pw2wannier90.x : pw2wannier90.o $(PPOBJS) $(MODULES) $(LIBOBJS)
diff --git a/PW/Makefile b/PW/Makefile
index a611bec..baddd6a 100644
--- a/PW/Makefile
+++ b/PW/Makefile
@@ -145,6 +145,7 @@ print_ks_energies.o \
 punch.o \
 pw_restart.o \
 pwcom.o \
+pw2blip.o \
 qvan2.o \
 rcgdiagg.o \
 rdiagh.o \
@@ -229,7 +230,8 @@ wannier_init.o \
 wannier_check.o \
 wannier_clean.o \
 wannier_occ.o \
-wannier_enrg.o 
+wannier_enrg.o \
+write_casino_wfn.o
 
 EEOBJS=../EE/libee.a
 QEMODS=../Modules/libqemod.a
diff --git a/PW/pwscf.f90 b/PW/pwscf.f90
index 7fbd607..f72582e 100644
--- a/PW/pwscf.f90
+++ b/PW/pwscf.f90
@@ -29,6 +29,8 @@ PROGRAM pwscf
 #endif
   !
   IMPLICIT NONE
+  CHARACTER(4) postfix
+  INTEGER itercount
   !
   !
 #ifdef __PARA
@@ -84,8 +86,11 @@ PROGRAM pwscf
      !
      CALL init_run()
      !
+     itercount=0
+     !
      main_loop: DO
         !
+        itercount=itercount+1
         !
         ! ... electronic self-consistentcy
         !
@@ -93,6 +98,16 @@ PROGRAM pwscf
         !
         IF ( .NOT. conv_elec ) CALL stop_run( conv_elec )
         !
+        write(postfix,'(i4.4)')itercount
+        CALL write_casino_wfn( &
+            .false., & ! gather
+            .true.,  & ! blip
+            1.0d0,   & ! multiplicity
+            .true.,  & ! binwrite
+            .true.,  & ! single_precision_blips
+            0,       & ! n_points_for_test
+            '.'//postfix)   ! postfix
+        !
         ! ... if requested ions are moved
         !
         CALL ions()
EOF

patch -p 1 --dry-run -f -p 1 < $PATCHFILE

if [ "$?" != "0" ] ; then
    echo
    echo "ERROR: Patch did not apply cleanly."
    echo "--> Very likely, the distribution has been changed without updating the script PP/$SCRIPTNAME"
    echo "--> Either apply the patch $PATCHFILE manually (it is fairly simple) or contact the author for help:"
    echo "-->     <Norbert@Nemec-online.de>"
    exit 1
fi

patch -p 1 -s -b < $PATCHFILE
rm -f $PATCHFILE

mv PP/pw2blip.f90 PW/pw2blip.f90
mv PP/pw2casino_write.f90 PW/write_casino_wfn.f90
install/makedeps.sh
