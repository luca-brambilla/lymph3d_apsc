#!/bin/bash

LYMPH_DIR="$(pwd)"

# METIS
cd metis-5.1.0
make config prefix=./INSTALLED
make
make install

# FMETIS
cd ../fmetis
mkdir -p build
cd build
rm -rf *
cmake .. -DMETIS_LIB="$(LYMPH_DIR)/metis-5.1.0/build/Linux-x86_64/INSTALLED/lib/libmetis.a" -DREAL=64 -DCMAKE_INSTALL_PREFIX=./INSTALLED
make install
ctest --rerun-failed --output-on-failure
