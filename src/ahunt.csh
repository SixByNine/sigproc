###############################################################################
#
# ahunt - script to dedisperse and search data using seek over a given DM range
#
# Also does accelerated search at each DM. MJK 2005
# 
# The script requires a filterbank file with the raw data and a "dmlist" file
# which is an ASCII list of the trial DMs to use in the search.
#
###############################################################################
set dedisperse = $bin/dedisperse
set seek       = $bin/seek
set accn       = $bin/accn
set source = $1
if ($source == "") then
	echo "usage: ahunt filestem (option)"
	exit
endif
if (! -e $source.fil) then
	echo "Filterbank file: $source.fil not found..."	
	exit
endif
if (! -e dmlist) then
	echo "File dmlist with DM ranges not found..."
	exit
endif
###############################################################################
echo "Hunting file $source.fil for periodicities... check back here later!"
echo "DM:" | awk '{printf "%s ", $1}'
set append = ""
foreach dm (`cat dmlist`)
	echo $dm | awk '{printf "%s ", $1}'
	$dedisperse $source.fil -d $dm > $source.tim
	$accn $source.tim $append > /dev/null
	set append = "-A"
end
echo DONE
exit
###############################################################################
