#!/bin/bash
#
# Shell for creating a tarball with DAMQT package for installation with cmake
#
# To run this shell appropriately, in the home directory of DAMQT execute
#
#    utils/prepare_tarball.sh
#
# and it will create a tarball in the upper directory 
# 
# Author: Rafael Lopez (rafael.lopez@uam.es)
# Date: March 2016
# =============================================================================
export directorio=${PWD##*/} 
# echo "directorio = " $directorio
if [ "$1" != "" ]; then
   filename=$1  
else
   fecha=$(date "+%Y%m%d")
#    echo $fecha
   filename=${directorio}_${fecha}
fi
# echo "filename = " $filename

echo " "
echo "-------------------------------------------------------------------------------------------------------------------------"
echo "Creates a tarball of " $directorio " in "$(dirname ${PWD})/$filename".tgz"
echo "-------------------------------------------------------------------------------------------------------------------------"
echo " "

tar --wildcards -cvzh --exclude="bin" --exclude="bin/*" --exclude="release" \
--exclude="release/*" --exclude="Release" --exclude="Release/*" --exclude="debug" --exclude="debug/*"  --exclude="Debug" \
--exclude="Debug/*" \
-f "$(dirname ${PWD})/$filename".tgz ../$directorio/*
