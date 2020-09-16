Install

# HDF 5.10.4

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz

./configure --prefix=/home/buildbot/local/hdf5-104-gcc102 --enable-fortran

make -j 4
make install
