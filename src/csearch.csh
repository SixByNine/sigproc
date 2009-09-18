###############################################################################
#
# csearch - script to do a companion search on a dedispersed time series
#
# The script requires a .tim file with the dedispersed data and an "aclist" file
# which is an ASCII list of the trial ACs and ADOTs to use in the search.
#
###############################################################################
set seek       = $bin/seek
set header     = $bin/header
set polyco     = $bin/polyco
set depolyco   = $bin/depolyco
set source     = $1
set pulsar     = $2
if ($source == "") then
	echo "usage: csearch filestem parstem"
	exit
endif
if (! -e $source.tim) then
	echo "Time series file: $source.tim not found..."	
	exit
endif
if (! -e a1list) then
	echo "File a1list with A1 ranges not found..."
	exit
endif
if (! -e $pulsar.par) then
	echo "TEMPO par file not found..."
	exit
endif
echo "Hunting file $source.tim for periodicities... check back here later"
echo "A1:" | awk '{printf "%s ", $1}'
set append = ""
set mjd = `$header $source.tim -tstart`
set tel = `$header $source.tim -telescope`
set frq = `$header $source.tim -fsamp`
foreach a1 (`cat a1list`)
	echo $a1 | awk '{printf "%s ", $1}'
	cp $pulsar.par $TEMPO/tzpar
	echo "A1 $a1"  >> $TEMPO/tzpar/$pulsar.par
	echo "F0 $frq" >> $TEMPO/tzpar/$pulsar.par
	$polyco $pulsar -mjd $mjd -site $tel -maxha 12 -nspan 60 > /dev/null
	$depolyco $source.tim polyco.dat > $pulsar.tim
	$seek $pulsar.tim -q $append -D$a1 -c4
	set append = "-A"
end
exit
