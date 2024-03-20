#!/bin/bash
git clone https://github.com/pdfpc/pdfpc.git
cd pdfpc
mkdir build

cd build/

cmake -DCMAKE_INSTALL_PREFIX="$(pwd)" ..
make
make install
