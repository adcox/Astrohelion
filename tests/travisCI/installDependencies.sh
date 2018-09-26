#!/bin/sh

# Install dependencies for Astrohelion on a Linux 64-bit machine
# All included files, libraries, and executables will be placed in
# the directory specified by INSTALL_DIR

# Save some (very lengthy) output to this file [doesn't seem to work...]
LOG="installLog.log"

# If you don't have permissions for this directory, run script with sudo
INSTALL_DIR="/usr/local"

# Settings, body data live in this folder
CONFIG_DIR="$HOME/.config/astrohelion"

# SPICE data lives in this directory
SPICE_DIR="$HOME/.config/astrohelion/data/data_SPICE"
SPK="de430.bsp"		# Planetary body ephemeris SPICE kernel
TLS="naif0012.tls"	# Time (leap second) SPICE kernel

# Figure out the OS
OS="unknown"
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
	OS="linux"
elif [[ "$unamestr" == "Darwin" ]]; then
	OS="darwin"
fi

mkdir -p deps
cd deps

## Boost
wget https://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz
tar -xzf boost*.tar.gz >> $LOG
cd boost*
echo "Installing BOOST"
./bootstrap.sh --prefix=$INSTALL_DIR
if [[ "$OS" == "darwin" ]]; then
	./b2 toolset=gcc-6 -d0 install --with-test
else
	./b2 -d0 install --with-test
fi
cd ..
rm boost*.tar.gz

## CSpice
if [[ "$OS" == "darwin" ]]; then
	wget http://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z
else 
	wget http://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
fi
tar -xzf cspice.tar.Z >> $LOG
cd cspice
echo "Installing CSPICE"
mkdir -p $INSTALL_DIR/include/cspice
mkdir -p $INSTALL_DIR/bin/cspice
cp -R include/* $INSTALL_DIR/include/cspice/
mv lib/cspice.a $INSTALL_DIR/lib/libcspice.a
mv lib/csupport.a $INSTALL_DIR/lib/libcsupport.a
mv exe/* $INSTALL_DIR/bin/cspice/
cd ..
rm cspice.tar.Z

## MatIO
wget https://sourceforge.net/projects/matio/files/matio/1.5.11/matio-1.5.11.tar.gz
tar -xzf matio-* >> $LOG
cd matio-*
echo "Installing MATIO"
./configure --prefix=$INSTALL_DIR
make >> $LOG
make install >> $LOG
cd ..
rm matio-1.5.11.tar.gz

## GSL
wget http://mirror.nexcess.net/gnu/gsl/gsl-2.4.tar.gz
tar -xzf gsl-*.tar.gz >> $LOG
cd gsl-*
echo "Installing GSL"
./configure --prefix=$INSTALL_DIR
make >> $LOG
make install >> $LOG
cd ..
rm gsl-*.tar.gz

if [[ "$OS" == "linux" ]]; then
	libtool --finish $INSTALL_DIR/lib/
fi

## Download and move SPICE data
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/$SPK
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/$TLS
mkdir -p $SPICE_DIR
mv $SPK $SPICE_DIR/
mv $TLS $SPICE_DIR/

## Eigen - Header only
wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
tar -xzf 3.3.3.tar.gz
mv eigen* eigen
mv eigen/Eigen $INSTALL_DIR/include/Eigen
mv eigen/unsupported $INSTALL_DIR/include/Eigen/unsupported
rm 3.3.3.tar.gz

## Configuration Files
cd ..
cp ../../misc/body_data.xml $CONFIG_DIR/body_data.xml
SETTINGS="$CONFIG_DIR/user_settings.xml"
echo "<?xml version=\"1.0\" encoding=\"utf-8\"?>" > $SETTINGS
echo "<astrohelion>" >> $SETTINGS
echo "  <spice>" >> $SETTINGS
echo "    <data_filepath>$SPICE_DIR/</data_filepath>" >> $SETTINGS # End with <filesep>
echo "    <time_kernel>$TLS</time_kernel>" >> $SETTINGS
echo "    <spk_kernel>$SPK</spk_kernel>" >> $SETTINGS
echo "  </spice>" >> $SETTINGS
echo "</astrohelion>" >> $SETTINGS

cp $SETTINGS $CONFIG_DIR/default_settings.xml

## Extras for the Travis-CI system
echo "Adding $INSTALL_DIR/lib to LD_LIBRARY_PATH"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/lib
LD_RUN_PATH=$LD_RUN_PATH:$INSTALL_DIR/lib

# Update the dynamic linker
sudo ldconfig