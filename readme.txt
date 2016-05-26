First compile sdsl-lite with:

git submodule init
git submodule update
cd sdsl-lite
sh ./install.sh

Then compile the rest with make at the project root.

Library files will be compiled into ./lib and
headers will be installed to ./include

To use the library, link against -lbdbwt -ldbwt -ldivsufsort64 -lsdsl
For example: 

g++ main.cpp -std=c++11 -L lib -I include -lbdbwt -ldbwt -ldivsufsort64 -lsdsl

The header include/BD_BWT_index.hh contains some documentation.
