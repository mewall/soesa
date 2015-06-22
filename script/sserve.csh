#!/bin/csh -f
setenv r <place soesa root directory here>
setenv b {$r}/bin
setenv d {$r}/dat
nohup $b/soesa -aacfg $d/aa.cfg -hashfile $d/pdf.hash -datafile $d/pdf.data \
         -servermode -verbose -fifo soesa.fifo -pdb user.pdb -out user.dat \
	 >& soesa.log &
