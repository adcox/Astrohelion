#!/bin/sh

# INSTALL_DIR is the directory where libraries, headers, and other resources are
# installed. You should choose a directory you have rwx priveleges for.
INSTALL_DIR="/opt/local"

# CONFIG_DIR is a directory where Astrohelion configuration files are saved. By
# default, the location is within the .config directory found in all Unix home 
# directory.
CONFIG_DIR="$HOME/.config/astrohelion"

# SPICE_DIR specifies the directory for SPICE data files. By default, this is 
# located in the same place as the configuration files. However, if you have 
# another directory with SPICE files, feel free to change this variable; the 
# settings need not be collocated with the SPICE files.
SPICE_DIR="$HOME/.config/astrohelion/spice"
SPK="de430.bsp"		# Planetary body ephemeris SPICE kernel
TLS="naif0012.tls"	# Time (leap second) SPICE kernel


# ##############################################################################
# 	DO NOT EDIT BELOW THIS LINE
# ##############################################################################

# Figure out the OS
OS="unknown"
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
	OS="linux"
elif [[ "$unamestr" == "Darwin" ]]; then
	OS="darwin"
fi

# Create directory for dependencies (may be deleted after installation)
mkdir -p deps
cd deps

## Boost
wget https://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz
tar -xzf boost*.tar.gz
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
tar -xzf cspice.tar.Z
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
tar -xzf matio-*
cd matio-*
echo "Installing MATIO"
./configure --prefix=$INSTALL_DIR
make
make install
cd ..
rm matio-1.5.11.tar.gz

## GSL
wget http://mirror.nexcess.net/gnu/gsl/gsl-2.4.tar.gz
tar -xzf gsl-*.tar.gz
cd gsl-*
echo "Installing GSL"
./configure --prefix=$INSTALL_DIR
make
make install
cd ..
rm gsl-*.tar.gz

if [[ "$OS" == "linux" ]]; then
	libtool --finish $INSTALL_DIR/lib/
fi

## Download and move SPICE data
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/$SPK
wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/$TLS
mkdir -p $SPICE_DIR
echo "Copying SPICE files to $SPICE_DIR"
mv $SPK $SPICE_DIR/
mv $TLS $SPICE_DIR/

## Eigen - Header only
wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
tar -xzf 3.3.3.tar.gz
mv eigen* eigen
mv eigen/Eigen $INSTALL_DIR/include/Eigen
mv eigen/unsupported $INSTALL_DIR/include/Eigen/unsupported
rm 3.3.3.tar.gz

## Body data XML file
mkdir -p $CONFIG_DIR
echo "Copying body data to $CONFIG_DIR"
cp misc/body_data.xml $CONFIG_DIR/body_data.xml

## Write settings to XML file
SETTINGS="$CONFIG_DIR/user_settings.xml"
echo "<?xml version=\"1.0\" encoding=\"utf-8\"?>" > $SETTINGS
echo "<astrohelion>" >> $SETTINGS
echo "  <spice>" >> $SETTINGS
echo "    <data_filepath>$SPICE_DIR/</data_filepath>" >> $SETTINGS # End with <filesep>
echo "    <time_kernel>$TLS</time_kernel>" >> $SETTINGS
echo "    <spk_kernel>$SPK</spk_kernel>" >> $SETTINGS
echo "  </spice>" >> $SETTINGS
echo "</astrohelion>" >> $SETTINGS

## Copy the settings file to the defaults
echo "Saving settings to $CONFIG_DIR"
cp $SETTINGS $CONFIG_DIR/default_settings.xml