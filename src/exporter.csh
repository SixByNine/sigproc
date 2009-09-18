#!/bin/csh
# exporter.csh - script to produce a tar file for exporting to other machines
set d = `pwd | awk -F/ '{print $NF}'`
set v = `awk '{print $1}' version.history | tail -2 | head -1`
set date = `date +"%A %B %e, %Y"`
echo "% THIS FILE SHOULD NOT BE EDITED -- EDIT documentation.tex" >! sigproc.tex
cat documentation.tex|sed s/"X.X"/"$v"/|sed s/"RELEASE"/"$date"/>>sigproc.tex
echo "#define SIGPROC_VERSION $v" > version.h
echo "      character*60 version" >! vers.inc
echo "      parameter(version='is part of SIGPROC version: $v')" >> vers.inc
cd ..
if ($d != "sigproc-$v") mv $d sigproc-$v
set d = sigproc-$v
if (-e $d.tar) rm -f $d.tar
if (-e $d.tar.gz) rm -f $d.tar.gz
tar cf $d.tar $d/version.history $d/makefile $d/configure $d/*.tex \
       $d/*.ps $d/*.c $d/*.C $d/*.h $d/*.l $d/*.y $d/*.csh $d/*.tcl $d/*.tk \
       $d/*.f $d/*.inc $d/*.pdf
gzip   $d.tar
echo "export version of $d ready..."
ls -l  $d.tar.gz
exit
