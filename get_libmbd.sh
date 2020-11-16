#!/bin/bash
git clone --depth 1 --branch 0.10.3  https://github.com/jhrmnn/libmbd.git
sed -i '29i#GSz: QE Support\nIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)\n  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Intstalling to local directory" FORCE)\nENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)' libmbd/CMakeLists.txt
cd ./libmbd/
cmake .
make install
if [[ -d lib64 ]]; then
  ln -s lib64 lib
fi
cd -
