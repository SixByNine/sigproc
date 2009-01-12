###############################################################################
#
# foldsignals - script to dedisperse and fold signals obtained from SEEK
#
# The script first creates a list of signals by running "best" (with user
# specified s/n threshold). Then, for each one, it attempts to optimize the
# search period by folding over different trial periods within +/- one 
# fourier bin of the search period. The best period is judged to be the 
# resulting profile with the biggest "bump" which is measured by its rms.
#
###############################################################################
set dedisperse = $bin/dedisperse
set fold       = $bin/fold
set grey       = $bin/grey
set decimate   = $bin/decimate
set header     = $bin/header
set quickplot  = $bin/quickplot
set reader     = $bin/reader
set best       = $bin/best
set ntrials    = 64
set nbands     = 8
set nsubints   = 8
set snmin      = 8
set nfold      = 999
set userbins   = 0
set usercand   = 0
set listonly   = 0
###############################################################################
set nargs = `echo $argv | wc -w`
if ($nargs == 0) then
    echo ""
    echo "foldsignals - script to dedisperse and fold signals from SEEK"
    echo ""
    echo "usage: foldsignals stem -{options}"
    echo ""
    echo "options:"
    echo ""
    echo "-nbands n      - number of sub-bands in output plot (def=8)"
    echo "-ntrials n     - number of trial periods to try (def=64)"
    echo "-nsubints n    - number of sub-integrations in output plot (def=8)"
    echo "-snmin s       - miniumum signal-to-noise ratio (def=8)"
    echo "-nbins n       - number of bins in output profile (def=P/tsamp)"
    echo "-ncand n       - process only candidate number n (def=all)"
    echo "-listonly      - just list the candidates from best"
    echo ""
    exit
endif
if ($nargs > 1) then
    set i = 2
    while ($i <= $nargs) 
        if ("$argv[$i]" == "-listonly") then
           set listonly = 1
	endif
        if ("$argv[$i]" == "-nbands") then
	   @ j = $i + 1
           set nbands = $argv[$j]
	endif
        if ("$argv[$i]" == "-ntrials") then
	   @ j = $i + 1
           set ntrials = $argv[$j]
	endif
        if ("$argv[$i]" == "-ncand") then
	   @ j = $i + 1
           set usercand = $argv[$j]
	endif
        if ("$argv[$i]" == "-nbins") then
	   @ j = $i + 1
           set userbins = $argv[$j]
	endif
	if ("$argv[$i]" == "-nsubints") then
	   @ j = $i + 1
           set nsubints = $argv[$j]
	endif
	if ("$argv[$i]" == "-snmin") then
	   @ j = $i + 1
           set snmin = $argv[$j]
	endif
	if ("$argv[$i]" == "-nfold") then
	   @ j = $i + 1
           set nfold = $argv[$j]
	   set snmin = 5.0
	endif
	@ i++
    end
endif
set stem = $argv[1]
if (! -e $stem.prd) then
    echo "$1 not found...."
    exit
endif
if (! -e $stem.fil) then
    echo "filterbank file $stem.fil not found...."
    exit
endif
###############################################################################
set cal = `$header $stem.fil -date`
set obstime = `header $stem.fil -tobs`
set ra = `$header $stem.fil -src_raj`
set de = `$header $stem.fil -src_dej`
set coords = "$ra$de"
set nchans = `$header $stem.fil -nchans`
set fch1 = `$header $stem.fil -fch1`
set foff = `$header $stem.fil -foff`
set fmid = `echo "$fch1+$foff*($nchans/2)" | bc -l`
set chan = `echo $foff | awk '{if ($1<0) {print -1.0*$1} else {print $1}}'`
set tsamp = `$header $stem.fil -tsamp`  
###############################################################################
if (-e $stem.1.ps.gz) rm -f $stem*.ps.gz
if (-e $stem.001.fld) rm -f $stem*.fld
$best $stem.prd -p -F$fmid -C$chan -s$snmin -z > /dev/null
set ncands = `cat $stem.lis | wc -l`
if ($nfold < 999) then
  echo "Filterbank file: $stem.fil Top $nfold signals from the search:"
else
  echo "Filterbank file: $stem.fil $ncands signals found with S/N > $snmin..."
endif    
echo "------------------------------------------------------------"
echo " Period (ms)    S/N     DM   DMID Nhit Fold P/Ptop    Ptop/P"
echo "------------------------------------------------------------"
cat $stem.lis | head -$nfold
echo "------------------------------------------------------------"
if ($listonly == 1) exit
set cand = 0
echo "Folding signal:" | awk '{printf "%s ",$0}'
foreach period (`awk '{print $1}' $stem.lis`)
  set dmidx = `grep $period $stem.lis | awk '{print $4}'`
  set dmval = `grep $period $stem.lis | awk '{print $3}'`
  set snr   = `grep $period $stem.lis | awk '{print $2}'`
  @ cand++
  if (($usercand > 0) && ($cand != $usercand)) goto next
  if ($cand > $nfold) goto done
  set srch  = sus$cand.sum
  echo $cand | awk '{printf "%s ",$0}'
  if (! -e $dmidx.sub) then
        $dedisperse $stem.fil -d $dmval -f $fch1 -b $nbands > $dmidx.sub
        $dedisperse $dmidx.sub -d 0 > $dmidx.tim
  endif
  set candfile = `echo "$stem $cand" | awk '{printf "%s.%03d.fld",$1,$2}'`
  echo "--------------------------------------------------------------" \
    >> $candfile
  echo "Signal number: $cand Period: $period ms DM: $dmval S/N: $snr" \
    >> $candfile
  if ($userbins <= 0) then
  set nbins = `echo "(1000.0*$period/$tsamp)"| bc -l | awk '{if ($1>128.0) print 128; if ($1<128.0) print int($1)}'`
  else
  set nbins = $userbins
  endif
  set fsearch = `echo 1000.0/$period | bc -l`
  set delta_f = `echo 1.0/$obstime | bc -l`
  set pfldmin = `echo "1000.0/($fsearch+$delta_f)" | bc -l`
  set pfldmax = `echo "1000.0/($fsearch-$delta_f)" | bc -l`
  if ($ntrials != 0) then
    set delta_p = `echo "($pfldmax-$pfldmin)/$ntrials.0" | bc -l`
  endif
  set p_trial = $pfldmin
  set bestrms = 0.0
  set bestprd = $period
  set i = 0
  while ($i < $ntrials)
    $fold $dmidx.tim -p $p_trial -n $nbins | tail -$nbins \
	  | awk '{print $2}' >! test.profile
    set rms     = `cat test.profile|awk '{ssq+=$1*$1;print sqrt(ssq)}'|tail -1`
    set result  = `echo "$rms $bestrms" | awk '{if ($1>$2) print "YES"}'`
    if ($result == "YES") then
	set bestrms = $rms
	set bestprd = $p_trial
    endif
    set p_trial = `echo "($p_trial+$delta_p)" | bc -l`
    @ i++
  end
  rm -f test.profile
  set period =  $bestprd
  # now fold the data using the best period and produce the plot
  set tsub = `echo "$obstime/$nsubints"| bc -l | awk '{printf "%.5f", $1}'`
  $fold $dmidx.tim -p $period -d $tsub -n $nbins -epn >subints.epn
  $grey subints.epn | awk '{print substr($0,16)}' >> $candfile
  $fold $dmidx.sub -p $period -n $nbins -epn > subbands.epn
  $decimate $dmidx.tim -T 1024 | $reader >! timeseries
  echo "SEEK search output for" >! asciiheader
  set tel = `$header $stem.fil -telescope`
  set mac = `$header $stem.fil -machine`
  set frame = `$header $stem.fil -frame`
  set nbits = `$header $stem.fil -nbits`
  echo "$tel $mac $frame data" >> asciiheader
  echo $stem.$cand >> asciiheader
  echo $coords >> asciiheader
  $header $stem.fil -tstart | awk '{printf "%.5f\n", $1}' >> asciiheader
  echo $cal  >> asciiheader
  echo "$fch1+($nchans*$foff)/2."|bc -l|awk '{printf "%.1f\n",$1}'>>asciiheader
  echo "($nchans*$foff)"         |bc -l|awk '{printf "%.1f\n",$1}'>>asciiheader
  echo $fch1 >> asciiheader
  echo $foff >> asciiheader
  echo $nchans >> asciiheader
  echo $nbits >> asciiheader
  echo $tsamp  >> asciiheader
  echo $obstime >> asciiheader
  echo $period >> asciiheader
  echo $dmval >> asciiheader
  set psf = `$quickplot $srch greyscale`
  rm -f $srch
  next:
end
###############################################################################
done:
echo "DONE!"
echo "--------------------------------------------------------------" >! $stem.txt
echo " Period (ms)    S/N     DM   DMID Nhit Fold P/Ptop    Ptop/P" >>$stem.txt
echo "--------------------------------------------------------------" >>$stem.txt
cat $stem.lis | head -$nfold                                        >>$stem.txt
cat $stem.*.fld                                             >> $stem.txt
echo "--------------------------------------------------------------" >>$stem.txt
###############################################################################
cleanup:
rm -f *.top *.frq *.monitor *.tim *.sub *.log *.bst *.sum grey *FFT* *FFA*
rm -f $stem sus* 
rm -f *core* *.fld *.lis 
rm -f asciiheader sub*.epn timeseries
rm -f $stem
gzip -f *.ps
echo "ASCII output in file: $stem.txt"
echo "The corresponding plot files are:"
ls $stem*.ps.gz
exit
###############################################################################
