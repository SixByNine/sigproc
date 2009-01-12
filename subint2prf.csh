#!/bin/csh
if ($1 == "") then
	echo "usage: subint2prf file"
	exit
endif
set mjd = `head -1 $1 | awk '{printf "%d",$2}'`
set f = $1
set n = `head -1 $f | awk '{print $8}'`
@ h = $n + 1
@ b = $n + 1
set n = `wc -l $f | awk '{print $1}'`
set i = 1
while ($h < $n)
head -$h $f | tail -$b > $mjd.$i.prf
@ i++
@ h = $h + $b
end
