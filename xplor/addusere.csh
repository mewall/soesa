#!/bin/csh -x
# 
# Script to link 0usere into XPLOR, enabling SOESA to be used via
#   the USER potential.
# 
# XPLOR environment variables need to be defined (new-sgi-ulogin.com)
# Object (.o) files from the installation must still be in $OBJ.
#
# First copy the object file 0usere-[platform]-[type].o to $OBJ dir.
# Then ensure that this script is executable via chmod +x addusere.csh
# Finally execute this script from a csh shell via ./addusere.csh
#
cc -c -o $OBJ/calcuser.o $XPLOR_MIPS calcuser.c
$COMS/link.com 
