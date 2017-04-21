#!/bin/sh

LOG="installLog.log"
INSTALL_DIR="/opt/local"

# Figure out the OS
OS="unknown"
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
	OS="linux"
elif [[ "$unamestr" == "Darwin" ]]; then
	OS="darwin"
fi

# Install Dependencies for Travis CI system
mkdir -p deps
cd deps

sudo apt-get -qq update

## Boost
wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz
tar -xvzf boost_1_62_0.tar.gz >> $LOG
cd boost_1_62_0
echo "Installing BOOST"
./bootstrap.sh --prefix=$INSTALL_DIR >> $LOG
if [[ "$OS" == "darwin" ]]; then
	sudo ./b2 toolset=gcc-6 -d0 install --with-filesystem --with-system --with-test >> $LOG
else
	sudo ./b2 -d0 install --with-filesystem --with-system --with-test >> $LOG
fi
cd ..

## CSpice
if [[ "$OS" == "darwin" ]]; then
	wget http://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z
else 
	wget http://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
fi
tar -xvzf cspice.tar.Z >> $LOG
cd cspice
echo "Installing CSPICE"
sudo mkdir -p $INSTALL_DIR/include/cspice
sudo mkdir -p $INSTALL_DIR/bin/cspice
sudo cp -R include/* $INSTALL_DIR/include/cspice/
sudo mv lib/cspice.a $INSTALL_DIR/lib/libcspice.a
sudo mv lib/csupport.a $INSTALL_DIR/lib/libcsupport.a
sudo mv exe/* $INSTALL_DIR/bin/cspice/
cd ..

## MatIO
wget https://sourceforge.net/projects/matio/files/matio-1.5.10.tar.gz
tar -xvzf matio-* >> $LOG
cd matio-*
echo "Installing MATIO"
./configure --prefix=$INSTALL_DIR >> $LOG
make >> $LOG
sudo make install >> $LOG
cd ..

## GSL
wget http://mirror.nexcess.net/gnu/gsl/gsl-2.2.1.tar.gz
tar -xvzf gsl-2.2.1.tar.gz >> $LOG
cd gsl-2.2.1
echo "Installing GSL"
./configure --prefix=$INSTALL_DIR >> $LOG
make >> $LOG
sudo make install >> $LOG
cd ..

echo "Adding /usr/local/lib to LD_LIBRARY_PATH"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
LD_RUN_PATH=$LD_RUN_PATH:/usr/local/lib

# Update the dynamic linker
sudo ldconfig

## Download and move SPICE data
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
mkdir -p ~/.config/astrohelion
mv de430.bsp data/data_SPICE/
mv naif0012.tls data/data_SPICE/
cp ../travis_settings.xml ~/.config/astrohelion/user_settings.xml
cp ../body_data.xml ~/.config/astrohelion/body_data.xml

## Eigen - Header only
wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
tar -xvzf 3.3.3.tar.gz
mv eigen* eigen
sudo mv eigen/Eigen $INSTALL_DIR/include/Eigen
sudo mv eigen/unsupported $INSTALL_DIR/include/Eigen/unsupported
