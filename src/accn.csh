###############################################################################
#
# accn - script to do an acceleration search on a dedispersed time series
#
# The script requires a .tim file with the dedispersed data and an "aclist" file
# which is an ASCII list of the trial ACs and ADOTs to use in the search.
#
###############################################################################
set seek       = $bin/seek
set source = $1
if ($source == "") then
	echo "usage: accn filestem (option)"
	exit
endif
if (! -e $source.tim) then
	echo "Time series file: $source.tim not found..."	
	exit
endif
if (! -e aclist) then
	echo "File aclist with AC ranges not found..."
	exit
endif
if (! -e adlist) then
	echo 0.0 > adlist
endif
echo "Hunting file $source.tim for periodicities... check back here later"
echo "AC/AD:" | awk '{printf "%s ", $1}'
set append = ""
foreach ac (`cat aclist`)
	echo $ac | awk '{printf "%s/", $1}'
	foreach ad (`cat adlist`)
		echo $ad | awk '{printf "%s ", $1}'
		$seek $source.tim -q -a$ac -d$ad $append $2
		set append = "-A"
	end
end
exit
