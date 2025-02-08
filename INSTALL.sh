#!/bin/bash

LYMPH_DIR="$(pwd)"

# METIS
cd metis-5.1.0
make config #prefix=.
make
sudo make install

# FMETIS
cd ../fmetis
mkdir -p build
cd build
rm -rf *
cmake .. -DMETIS_LIB="/usr/local/lib/libmetis.a" -DREAL=64 #-DCMAKE_INSTALL_PREFIX=./installed
sudo make install
ctest --rerun-failed --output-on-failure