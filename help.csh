#!/bin/csh
# updates help files for inclusion in LaTeX manual
set program = `echo $1 | awk -F/ '{print $NF}'`
echo "\begin{verbatim}" >! $program.help
echo "% $program help" >> $program.help
$1 help >> $program.help
# strip off the last (blank) line of the help file 
set nlines = `cat $program.help | wc -l`
@ nlines = $nlines - 1
head -$nlines $program.help >! $program.temp
mv $program.temp $program.help
echo "\end{verbatim}" >> $program.help
chmod 666 $program.help
