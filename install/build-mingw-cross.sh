#!/bin/sh
# automated build script to build windows installers from the lammps sources

MINGW_BUILD_DIR=${HOME}/mingw-cross
NUMCPU=${NUMCPU-1}
echo X-compiling Quantum ESPRESSO for Windows in ${MINGW_BUILD_DIR}

for d in "${PWD}" "${PWD%/install}" "$1"
do \
    if test -d "${d}/.git"
    then
        ESPRESSO_PATH="${d}"
        break
    fi
done

if test -z "${ESPRESSO_PATH}"
then
    echo "'${PWD}' is not a suitable working directory"
    exit 1
fi

# clean up leftovers from an old build and rebuild directories
for d in qe-{serial,mpich2}-{32,64} espresso-current qe-docs
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done
mkdir -p ${MINGW_BUILD_DIR}/qe4win

pushd "${ESPRESSO_PATH}"

git archive -v --format=tar --prefix=espresso-current/ HEAD \
    | tar -C ${MINGW_BUILD_DIR} -xvf -
popd

mkdir -p ${MINGW_BUILD_DIR}/qe-serial-32
mkdir -p ${MINGW_BUILD_DIR}/qe-serial-64
mkdir -p ${MINGW_BUILD_DIR}/qe-mpich2-32
mkdir -p ${MINGW_BUILD_DIR}/qe-mpich2-64

pushd ${MINGW_BUILD_DIR}/espresso-current

# build and collect various pieces of documentation
./configure
make doc
pushd Doc
htmldoc --batch qe-input-ref.book
cp -v developer_man.pdf ${MINGW_BUILD_DIR}/qe-docs/QE_DeveloperManual.pdf
cp -v user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/QE_UserGuide.pdf
cp -v plumed_quick_ref.pdf brillouin_zones.pdf constraints_HOWTO.pdf ${MINGW_BUILD_DIR}/qe-docs/
cp -v QE-logo.jpg qe-input-ref.html ${MINGW_BUILD_DIR}/qe-docs/
cp -v release-notes ${MINGW_BUILD_DIR}/qe-docs/Release-Notes.txt
popd
cp -v License ${MINGW_BUILD_DIR}/qe-docs/License.txt
cp -v README ${MINGW_BUILD_DIR}/qe-docs/README.txt
cp -v CPV/Doc/user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/CPV_UserGuide.pdf
cp -v NEB/Doc/user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/NEB_UserGuide.pdf
cp -v PW/Doc/user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/PW_UserGuide.pdf
cp -v PP/Doc/user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/PP_UserGuide.pdf
cp -v PHonon/Doc/user_guide.pdf ${MINGW_BUILD_DIR}/qe-docs/PHonon_UserGuide.pdf
cp -v PHonon/Doc/developer_man.pdf ${MINGW_BUILD_DIR}/qe-docs/PHonon_DeveloperManual.pdf
cp -v atomic/Doc/pseudo-gen.pdf PP/Doc/eps_man.pdf TDDFPT/Doc/turboTDDFT-CPC.pdf ${MINGW_BUILD_DIR}/qe-docs/
unix2dos ${MINGW_BUILD_DIR}/qe-docs/*.txt
unix2dos ${MINGW_BUILD_DIR}/qe-docs/*.html

make distclean

# serial 32-bit
mingw32-configure LIBDIRS=$PWD/install/blas-win32 FFLAGS="-O3 -march=core2" CFLAGS="-O3 -march=core2"
make -j${NUMCPU} all || exit 1
make w90 || exit 1

STRIP=$(eval `rpm --eval %{mingw32_env}`; echo $STRIP)
pushd bin
for s in *.x
do \
    cp -v $s ${MINGW_BUILD_DIR}/qe-serial-32/$s.exe
    ${STRIP} -g ${MINGW_BUILD_DIR}/qe-serial-32/$s.exe
done
popd
cp install/blas-win32/libopenblas.dll ${MINGW_BUILD_DIR}/qe-serial-32

make distclean

# serial 64-bit
mingw64-configure FFLAGS="-O3 -march=core2" CFLAGS="-O3 -march=core2" # LIBDIRS=$PWD/install/blas-win64
make -j${NUMCPU} all || exit 1
make w90 || exit 1

STRIP=$(eval `rpm --eval %{mingw64_env}`; echo $STRIP)
pushd bin
for s in *.x
do \
    cp -v $s ${MINGW_BUILD_DIR}/qe-serial-64/$s.exe
    ${STRIP} -g ${MINGW_BUILD_DIR}/qe-serial-64/$s.exe
done
popd
#cp install/blas-win64/libopenblas.dll ${MINGW_BUILD_DIR}/qe-serial-64

make distclean

# mpich2 32-bit
mingw32-configure LIBDIRS=$PWD/install/blas-win32 FFLAGS="-O3 -march=core2" CFLAGS="-O3 -march=core2" \
    MPI_LIBS="-L$PWD/install/mpich2-win32/lib -lfmpi -lmpi"
make -j${NUMCPU} MANUAL_DFLAGS="-I$PWD/install/mpich2-win32/include" all || exit 1
make w90 || exit 1

STRIP=$(eval `rpm --eval %{mingw32_env}`; echo $STRIP)
pushd bin
for s in *.x
do \
    cp -v $s ${MINGW_BUILD_DIR}/qe-mpich2-32/$s.exe
    ${STRIP} -g ${MINGW_BUILD_DIR}/qe-mpich2-32/$s.exe
done
popd
cp install/blas-win32/libopenblas.dll ${MINGW_BUILD_DIR}/qe-mpich2-32

make distclean

# mpich2 64-bit

mingw64-configure FFLAGS="-O3 -march=core2" CFLAGS="-O3 -march=core2" \
    MPI_LIBS="-L$PWD/install/mpich2-win64/lib -lfmpi -lmpi"
# LIBDIRS=$PWD/install/blas-win64 
make -j${NUMCPU} MANUAL_DFLAGS="-I$PWD/install/mpich2-win64/include" all || exit 1
make w90 || exit 1
make -C W90/doc/user_guide
make -C W90/doc/tutorial
cp -v W90/doc/user_guide.pdf W ${MINGW_BUILD_DIR}/qe-docs/W90_UserGuide.pdf
cp -v W90/doc/tutorial.pdf W ${MINGW_BUILD_DIR}/qe-docs/W90_Tutorial.pdf

STRIP=$(eval `rpm --eval %{mingw64_env}`; echo $STRIP)
pushd bin
for s in *.x
do \
    cp -v $s ${MINGW_BUILD_DIR}/qe-mpich2-64/$s.exe
    ${STRIP} -g ${MINGW_BUILD_DIR}/qe-mpich2-64/$s.exe
done
popd
#cp install/blas-win64/libopenblas.dll ${MINGW_BUILD_DIR}/qe-mpich2-64

make distclean

popd

pushd ${MINGW_BUILD_DIR}
TOOLDIR=espresso-current/install

# make this the real version later on
verstr=5.1.svn$(date +%Y%m%d)
cp ${TOOLDIR}/espresso.nsis .
cp ${TOOLDIR}/EnvVarUpdate.nsh .

# determine os vendor and release for installer tweaks.
vendor=$(grep  release /etc/issue | cut -d \  -f 1)
release=$(grep  release /etc/issue | cut -d \  -f 3)
arch=$(uname -m)

# build installers
LIBGCC=libgcc_s_sjlj-1.dll
makensis -DMINGW=/usr/i686-w64-mingw32/sys-root/mingw/bin/ -DBIT=32 \
    -DVARIANT=serial -DVERSION=${verstr} -DLIBGCC=${LIBGCC} espresso.nsis
makensis -DMINGW=/usr/i686-w64-mingw32/sys-root/mingw/bin/ -DBIT=32 \
    -DVARIANT=mpich2 -DVERSION=${verstr} -DLIBGCC=${LIBGCC} espresso.nsis

# Fedora 19 ships with GCC-4.8.x which has different exception handling
# on 64-bit windows and thus uses a different name for libgcc
if [ "$vendor" = "Fedora" ] && [ $release -ge 19 ]
then
    LIBGCC=libgcc_s_seh-1.dll
fi
makensis -DMINGW=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/ -DBIT=64 \
    -DVARIANT=serial -DVERSION=${verstr} -DLIBGCC=${LIBGCC} espresso.nsis
makensis -DMINGW=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/ -DBIT=64 \
    -DVARIANT=mpich2 -DVERSION=${verstr} -DLIBGCC=${LIBGCC} espresso.nsis

exit 0
