#!/bin/sh

LOG="installLog.log"

# Install Dependencies for Travis CI system
mkdir -p deps
cd deps

sudo apt-get -qq update

# Boost
wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz
tar -xvzf boost_1_62_0.tar.gz >> $LOG
cd boost_1_62_0
echo "Installing BOOST"
./bootstrap.sh >> $LOG
sudo ./b2 -d0 install --with-filesystem --with-system --with-test >> $LOG 
cd ..

# CSpice
wget http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
tar -xvzf cspice.tar.Z >> $LOG
cd cspice
echo "Installing CSPICE"
sudo mkdir /usr/local/include/cspice
sudo mkdir /usr/local/bin/cspice
sudo cp -R include/* /usr/local/include/cspice/
sudo mv lib/cspice.a /usr/local/lib/libcspice.a
sudo mv lib/csupport.a /usr/local/lib/libcsupport.a
sudo cp exe/* /usr/local/bin/cspice/
cd ..

# MatIO
wget https://sourceforge.net/projects/matio/files/matio/1.5.3/matio-1.5.3.tar.gz
tar -xvzf matio-1.5.3.tar.gz >> $LOG
cd matio-1.5.3
echo "Installing MATIO"
./configure >> $LOG
make >> $LOG
sudo make install >> $LOG
cd ..

# GSL
wget http://mirror.nexcess.net/gnu/gsl/gsl-2.2.1.tar.gz
tar -xvzf gsl-2.2.1.tar.gz >> $LOG
cd gsl-2.2.1
echo "Installing GSL"
./configure >> $LOG
make >> $LOG
sudo make install >> $LOG
cd ..

echo "Adding /usr/local/lib to LD_LIBRARY_PATH"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
LD_RUN_PATH=$LD_RUN_PATH:/usr/local/lib

# Update the dynamic linker
sudo ldconfig

# Download and install SPICE data
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
mkdir -p ~/.config/astrohelion
mv de430.bsp ~/.config/astrohelion/
mv naif0012.tls ~/.config/astrohelion/
cp ../travis_settings.xml ~/.config/astrohelion/user_settings.xml
cp ../body_data.xml ~/.config/astrohelion/body_data.xml
