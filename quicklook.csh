# QUICKLOOK : a csh script to process fast-sampled pulsar data 
if ($1 == "" || $1 == "help") then
    echo ""
    echo "usage: quicklook <filename> -{options}"
    echo ""
    echo "options:"
    echo ""
    echo "-read time   - read and process only time (s) of data"
    echo "-skip time   - skip the first time (s) of data"
    echo "-addc nchans - add nchans chans together before dedispersion (def=none)"
    echo "-addt nsamps - add nsamps samples together before dedispersion (def=none)"
    echo "-nsints   n  - specify number of time sub-integrations (def=8)"
    echo "-nbands   n  - specify number of frequency sub-bands (def=8)"
    echo "-period p    - fold data at constant period (ms) (def=polyco)"
    echo "-dm dmvalue  - dedisperse using dmvalue (pc cm-3) (def=polyco)"
    echo "-clip value  - clip samples that deviate more than value*rms (def=noclipping)"
    echo "-nbins n     - fold data using n bins (def=128)"
    echo "-left phase  - specify left-hand phase window (def=0.0)"
    echo "-right phase - specify right-hand phase window (def=0.0)"
    echo "-singlepulse - set number of subints to be the number of pulses"
    echo "-psr name    - fix source name for running TEMPO (def=filestem)"
    echo "-mypolyco    - use an existing polyco file (def=make one)"
    echo "-fakescale   - plot channels/bands using fake grey scale (def=real)"
    echo "-clean       - remove large files before and after (def=keep them)"
    echo "-wipe        - remove large files beforehand (def=keep them)"
    echo "-ext string  - add a file extension string to output files (def=none)"
    echo ""
    exit
endif

# set up some default options
set addt = 1
set addc = 1
set clip = ""
set readopt = ""
set skipopt = ""
set sints = 8
set singlepulse = 0
set bands = 8
set period = "polyco.dat"
set dm = -1.0
set binopt = "-n 128"
set lopt = ""
set ropt = ""
set ext = ""
if (-e phasestart) rm -f phasestart
set polyname = ""
set mypolyco = 0
set clean = 0
set wipe = 0
set plotopt = "greyscale"

# check out state of input file and get a polyco name
set filename = $1
if (-e $filename) then
	  set polyname  = `echo $filename | awk -F. '{print substr($1,2,7)}'`
else
	  echo "file: $filename not found..."
 	  exit
endif

# check on free diskspace and exit if less remaining than the current file

set freespace = `df -k . | tail -1 | awk '{print $4*1024.0}'`
set filespace = `ls -ln $filename | awk '{print $5*1.0}'`
set full = \
`echo "$filespace $freespace"  | awk '{if ($1>$2) {print 1} else {print 0}}'`
set full = 0
if ($full == "1") then
	echo "not enough free space..."
	df -k .
	exit
endif

# parse the rest of the command line
set i = 1
set skipping = 1
foreach item ($argv)
	@ next = $i + 1
	switch ($item)
	case -read:
		set readopt = "-r $argv[$next]"
		set skipping = 1
		breaksw
	case -skip:
		set skipopt = "-s $argv[$next]"
		set skipping = 1
		breaksw
	case -nsints:
		set sints = $argv[$next]
		set skipping = 1
		breaksw
	case -singlepulse:
		set singlepulse = 1
		set skipping = 0
		breaksw
	case -nbands:
		set bands = $argv[$next]
		set skipping = 1
		breaksw
	case -addc:
		set addc = $argv[$next]
		set skipping = 1
		breaksw
	case -addt:
		set addt = $argv[$next]
		set skipping = 1
		breaksw
	case -clip:
		set clip = "-c $argv[$next]"
		set skipping = 1
		breaksw
	case -period:
		set period = $argv[$next]
		set skipping = 1
		breaksw
	case -left:
		set lopt = "-l $argv[$next]"
                echo $argv[$next] >! phasestart
		set skipping = 1
		breaksw
	case -right:
		set ropt = "-r $argv[$next]"
		set skipping = 1
		breaksw
	case -dm:
		set dm = $argv[$next]
		set skipping = 1
		breaksw
	case -nbins:
		set binopt = "-n $argv[$next]"
                if ($argv[$next] == "opt") set binopt = ""
		set skipping = 1
		breaksw
	case -psr:
		set source = $argv[$next]
		set skipping = 1
		breaksw
	case -mypolyco:
		set mypolyco = 1
		set skipping = 0
		breaksw
	case -clean:
		set clean = 1
                set wipe = 1
		set skipping = 0
		breaksw
	case -wipe:
		set wipe = 1
		set skipping = 0
		breaksw
	case -fakescale:
		set plotopt = ""
		set skipping = 0
		breaksw
	case -ext:
		set ext = $argv[$next]
		set skipping = 1
		breaksw
	default:
		if ($skipping == 0) then
                  quicklook help
		  echo "ERROR: command-line argument $item not recognized"
		  exit
		endif
		set skipping = 0
		breaksw
	endsw
	set i = $next
end

# wipe output data files before starting if requested
if ($wipe == "1") then
  if ($mypolyco == 0) rm -f polyco.dat
  rm -f $filename.fil $filename.sub $filename.tim subbands.epn subints.epn
endif

	set mjd = `filterbank $filename | header -tstart`
# establish whether we need to generate a polyco.dat file...
if ( ($period == polyco.dat) && ($mypolyco == 0) ) then	
	# get mjd and try to generate polyco
	echo "polyco $polyname -mjd $mjd"
	polyco $polyname -mjd $mjd  >! polyco.log
	set psize = `ls -l polyco.dat | awk '{print $5}'`
	if ($psize == 0) then
		echo "ERROR creating polyco.dat...."
		cat polyco.log
		exit
	endif
endif

# get the dm from the polyco if not already set by user
if ($dm == -1.0) set dm = `head -1 polyco.dat | awk '{print $5}'`

# produce filterbank file if necessary
if (! -e $filename.fil) then
   if (($addc>1) || ($addt>1)) then
	echo "filterbank $filename -sumifs $readopt $skipopt | decimate -t $addt -c $addc -n 32 > $filename.fil"
	      filterbank $filename -sumifs $readopt $skipopt | decimate -t $addt -c $addc -n 32 > $filename.fil
   else 
	echo "filterbank $filename -sumifs $readopt $skipopt > $filename.fil"
        filterbank $filename -sumifs $readopt $skipopt > $filename.fil
   endif
else
	echo "filterbank file: $filename.fil found"
endif

# make a little ascii header for later reading by the plot program
header $filename.fil >! $filename.hdr

# check whether we need to dedisperse a different number of subbands
if (-e $filename.sub) then
	set test = `header $filename.sub -nchans`
	if ($test != $bands) rm -f $filename.sub
else
	rm -f $filename.sub
endif

# check to see if user has asked for all channels to be folded
if ($bands == "all") then
        if (-e $filename.sub) rm -f $filename.sub
        # just make a link to the filterbank file!
        ln -s $filename.fil $filename.sub
else
# form dedispersed sub-bands if necessary
if (! -e $filename.sub) then
	echo "dedisperse $filename.fil -d $dm -b $bands $clip> $filename.sub"
	      dedisperse $filename.fil -d $dm -b $bands $clip> $filename.sub
endif
endif

# fold subbands
echo "fold $filename.sub -p $period $binopt $lopt $ropt -epn >! subbands.epn"
      fold $filename.sub -p $period $binopt $lopt $ropt -epn >! subbands.epn

# produce time series file if necessary (n.b sub bands only need to be added)
if (! -e $filename.tim) then
	echo "dedisperse $filename.sub -d $dm > $filename.tim"
	      dedisperse $filename.sub -d $dm > $filename.tim
else
	echo "time series file: $filename.tim found"
endif

# produce 1024-sample average of time series
        echo "decimate $filename.tim -T 1024 | reader >! timeseries"
              decimate $filename.tim -T 1024 | reader >! timeseries

# grab relevant parameters from the filterbank data header
set telescope = `header $filename.fil -telescope`
set machine   = `header $filename.fil -machine`
set fch1      = `header $filename.fil -fch1`
set foff      = `header $filename.fil -foff`
set nchans    = `header $filename.fil -nchans`
set nbits     = `header $filename.fil -nbits`
set tsamp     = `header $filename.fil -tsamp`
set tobs      = `header $filename.fil -tobs`
set fcent     = `echo "$fch1 $foff $nchans" | awk '{print $1+$2*($3-1.0)/2.0}'`
set bandwidth = `echo "$foff $nchans" | awk '{print sqrt($1*$1)*$2}'`
set tdump     = `echo "$tobs $sints" | awk '{print $1/$2}'`
if ($singlepulse == 1) set tdump = 1
set cal       = `grep "Gregorian" $filename.hdr | awk '{print $NF}'`

# put these in asciiheader file for use by the plot program
cat << END    >! asciiheader
$telescope
$machine  
$filename 
$polyname 
$mjd      
$cal      
$fcent    
$bandwidth
$fch1     
$foff     
$nchans   
$tsamp    
$tobs     
$period   
$dm       
END

# fold time series
echo "fold $filename.tim -p $period $binopt -d $tdump $lopt $ropt -epn >! subints.epn"
      fold $filename.tim -p $period $binopt -d $tdump $lopt $ropt -epn >! subints.epn

# create postscript file summary and master epn file 
quickplot $plotopt
mv $filename.ps $filename$ext.ps
cat subbands.epn subints.epn > $filename$ext.epn
if (-e polyco.dat) cp polyco.dat $filename$ext.pco
set date = `date`
echo "all done... $date"
echo "plotepn $filename$ext.epn"
echo "gv $filename$ext.ps"
echo "more $filename$ext.pco"

# wipe output data files unless requested otherwise
if ($clean == "1") then
  if ($mypolyco == 0) rm -f polyco.dat
  rm -f $filename.fil $filename.sub $filename.tim subbands.epn subints.epn
endif

exit
