#!/bin/csh
# makes a dummy TEMPO par and tim file to for depolyco
if ($1 == "") then
    echo "usage: makedummy timfile"
    exit
endif

set parfile = "0000+0000.par"
echo "PSR 0000+0000" >! $parfile
header $1 -src_raj | awk '{print "RAJ "$1}'       >> $parfile
header $1 -src_dej | awk '{print "DECJ "$1}'      >> $parfile
header $1 -tsamp   | awk '{print "F0  "1.0e6/$1}' >> $parfile
header $1 -tstart  | awk '{print "PEPOCH "$1}'    >> $parfile
echo "DM 0.0"                                     >> $parfile
