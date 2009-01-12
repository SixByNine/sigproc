#!/bin/csh
# script to output which files are different
foreach file (`ls`)
    if (-e $1/$file) then
    set diff = `diff -q $file $1 | awk '{print $NF}'`
    if ($diff == "differ") then
	echo $file
    endif
    else
    echo $file does not exist in $1
    endif
end
exit
