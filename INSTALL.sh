#!/bin/bash

LYMPH_DIR="$(pwd)"

# METIS
cd metis-5.1.0
make config prefix=.
make install

# FMETIS
cd ../fmetis
mkdir -p build
#export FMETIS_DIR
cd build
cmake .. -DMETIS_LIB="$LYMPH_DIR/metis-5.1.0/build/Linux-x86_64/lib/libmetis.a" -DREAL=64 -DCMAKE_INSTALL_PREFIX=./installed
make install
ctest