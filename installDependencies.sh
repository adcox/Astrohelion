#!/bin/sh

# Install Dependencies
mkdir deps
cd deps

sudo apt-get -qq update
sudo apt-get install gsl matio

wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz
tar -xvzf boost_1_62_0.tar.gz
cd boost_1_62_0 && ./boostrap.sh && sudo ./b2 install
cd ..

wget http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
tar -xvzf cspice.tar.Z
cd cspice
cp -R include/* /usr/local/include/cspice/
mv lib/cspice.a /usr/local/lib/libcspice.a
mv lib/csupport.a /usr/local/lib/libcsupport.a
cp exe/* /usr/local/bin/cspice/
cd..