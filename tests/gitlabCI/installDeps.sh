#!/bin/sh

# Install dependencies for Astrohelion on a Linux 64-bit machine
# All included files, libraries, and executables will be placed in
# the directory specified by INSTALL_DIR

LOG="installLog.log"
INSTALL_DIR="/builds/adcox/Astrohelion/deps/installed"

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
mkdir -p $INSTALL_DIR

## Boost
if [ ! -f "$INSTALL_DIR/include/boost" ]; then
	wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz
	tar -xzf boost_1_62_0.tar.gz >> $LOG
	cd boost_1_62_0
	echo "Installing BOOST"
	./bootstrap.sh --prefix=$INSTALL_DIR >> $LOG
	if [[ "$OS" == "darwin" ]]; then
		./b2 toolset=gcc-6 -d0 install --with-test >> $LOG
	else
		./b2 -d0 install --with-test >> $LOG
	fi
	cd ..
	rm boost_1_62_0.tar.gz
fi

## CSpice
if [ ! -d "%INSTALL_DIR/include/cspice" ]; then
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
fi

## MatIO
if [ ! -e "INSTALL_DIR/include/matio.h" ]; then
	wget https://sourceforge.net/projects/matio/files/matio-1.5.10.tar.gz
	tar -xzf matio-* >> $LOG
	cd matio-*
	echo "Installing MATIO"
	./configure --prefix=$INSTALL_DIR >> $LOG
	make >> $LOG
	make install >> $LOG
	cd ..
	rm matio-1.5.10.tar.gz
fi

## GSL
if [ ! -d "$INSTALL_DIR/include/gsl" ]; then
	wget http://mirror.nexcess.net/gnu/gsl/gsl-2.2.1.tar.gz
	tar -xzf gsl-2.2.1.tar.gz >> $LOG
	cd gsl-2.2.1
	echo "Installing GSL"
	./configure --prefix=$INSTALL_DIR >> $LOG
	make >> $LOG
	make install >> $LOG
	cd ..
	rm gsl-2.2.1.tar.gz
fi

if [[ "$OS" == "linux" ]]; then
	libtool --finish $INSTALL_DIR/lib/
fi

## Download and move SPICE data
if [ ! -d /root/.config/astrohelion ]; then
	wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
	wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
	mkdir -p /root/.config/astrohelion
	mkdir -p /builds/adcox/Astrohelion/data/data_SPICE
	mv de430.bsp /builds/adcox/Astrohelion/data/data_SPICE
	mv naif0012.tls /builds/adcox/Astrohelion/data/data_SPICE
fi

## Eigen - Header only
if [ ! -d "$INSTALL_DIR/include/Eigen" ]; then
	wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
	tar -xzf 3.3.3.tar.gz
	mv eigen* eigen
	mv eigen/Eigen $INSTALL_DIR/include/Eigen
	mv eigen/unsupported $INSTALL_DIR/include/Eigen/unsupported
	rm 3.3.3.tar.gz
fi

## Configuration Files
cp ../settings.xml /root/.config/astrohelion/user_settings.xml
cp ../body_data.xml /root/.config/astrohelion/body_data.xml

## Final print out for debugging purposes
echo "ls -la /root/.config/astrohelion"
ls -la /root/.config/astrohelion
echo "ls -la $INSTALL_DIR"
ls -la $INSTALL_DIR
echo "ls -la $INSTALL_DIR/include"
ls -la $INSTALL_DIR/include
echo "ls -la $INSTALL_DIR/lib"
ls -la $INSTALL_DIR/lib

# Navigate back to the parent directory
cd $INSTALL_DIR
