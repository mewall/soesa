#!/bin/csh -f
setenv r <place soesa root directory here>
setenv b {$r}/bin
setenv d {$r}/dat
$b/soesa -aacfg $d/aa.cfg -hashfile $d/pdf.hash -datafile $d/pdf.data \
         -eval -pdb $1 -out $2
