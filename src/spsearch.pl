#!/usr/bin/perl
# 
# Runs a basic SP search for one beam. 
#
# Sarah B. Spolaor
# 20 Jul 2010
#
#
if ($#ARGV < 0){
    printhelp();
    exit;
}
#$calcdel   = "$gpipedir/calcdelay.pl";


#################
### (1) Set defaults and read input values.
#################
@file = ();
#$file[0] = "";
$lodm = 0;
$hidm = 500;
$sigthresh = 6;
$widthresh = 30;
$declim = "512";
$dmlim = 3;
$killmask = "";
$intwidth = 40; #microseconds
$intwidth = 1.25; # 25% smear tolerance
$i = 0;
while ($i <= $#ARGV){
    if ($ARGV[$i] eq '-h') {           # Help call
	printhelp();
	exit;
    } elsif (-e $ARGV[$i]) {           # File
	push(@file,$ARGV[$i]);
	$i++;
    } elsif ($ARGV[$i] eq '-mindm') {  # DM searching minimum
	$lodm = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-maxdm') {  # DM searching maximum
	$hidm = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-sig') {    # Threshhold for pulse detection
	$sigthresh = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-wid') {    # Threshhold for bin tolerance
	$widthresh = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-box') {    # Decimation limit
	$declim = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-zapdm') {    # RFI DM cutoff
	$dmlim = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-smear') {      # RFI DM cutoff
	$intwidth = $ARGV[++$i];
	$i++;
    } elsif ($ARGV[$i] eq '-k') {      # Channel mask
	$killmask = $ARGV[++$i];
	if (!-e $killmask){
	    print STDERR "\n\tERROR:\n\tKill mask $killmask does not exist.\n\n";
	    exit(1);
	}
	$i++;
    } else {
	print STDERR "\n\n\tERROR:\n\tCannot interpret input \"$ARGV[$i]\" on command line.\n\t\tSee -h menu for more information\n\n";
	exit(2);
    }
}
if ($file[0] eq '' || $#ARGV < 0){
    print STDERR "\n\tERROR:\n\tPlease supply a valid filterbank file on the command line.\n\n";
    exit(3);
}
if ($lodm > $hidm){
    print STDERR "\n\tERROR:\n\tMinimum DM ($lodm) > Maximum DM ($hidm)\n\n";
    exit(4);
}
if ($dmcut >= $hidm){
    print STDERR "\n\tERROR:\n\tRFI DM cutoff ($dmlim) >= Maximum DM ($hidm)\n\n";
    exit(5);
}    




# Report inputs and run search for each file
foreach $infile (@file){
    print STDERR "\n*********************************************\n";
    print STDERR "   Processing $infile\n";
    printf STDERR "  %12s  $lodm\n","Min DM:     ";
    printf STDERR "  %12s  $hidm\n","Max DM:     ";
    printf STDERR "  %12s  $sigthresh\n","SNR thresh: ";
    printf STDERR "  %12s  $widthresh samples\n","Width tol.: ";
    printf STDERR "  %12s  $declim samples\n","Max boxcar: ";
    printf STDERR "  %12s  $dmlim\n","DM cutoff:  ";
    print STDERR "*********************************************\n";

    @arr = split(/\//,$infile);
    $basename = $arr[-1];
    $origbasename = $basename;
    $basename =~ s/\.fil//;


#################
### (2) Run dedisperse_all
#################
    print STDERR "\nRunning dedisperse_all...\n";
    if ($killmask eq ''){
	system "dedisperse_all -d $lodm $hidm -g 1000000 -tol $intwidth $infile -G -wid $widthresh -sig $sigthresh -dec $declim -cut $dmlim -mb -l";
    } else { ###!!! REMOVE THE -i IN THESE TWO DED_ALL STATEMENTS
	system "dedisperse_all -d $lodm $hidm -g 1000000 -tol $intwidth $infile -G -wid $widthresh -sig $sigthresh -dec $declim -cut $dmlim -mb -l -k $killmask";
    }


#################
### (3) Check for candidates.
#################
    $ncands = `ls -1d *$basename*.pulse | wc -l`;
    chomp($ncands);
    if ($ncands == 0) {
	print STDERR " * * * No SP candidates in file. * * *";
	next;
    } else {
	print STDERR "\n ! ! ! Found $ncands candidates in file.\n ! ! ! Creating candidate plots...\n";
    }


#################
### (4) Read DM list.
#################
    @dmlist = ();
    $i=1;
    open(DMS, "$origbasename.dmlog");
    while (<DMS>){
	chomp;
	@temp = split(/\s+/,$_);
	push(@dmlist, "$temp[1]");
	$i++;
    }
    close(DMS);


#################
### (5) Set up plot links + make plots.
#################
    open(PLOTCANDS, "ls -1d *.pulse |");
    while(<PLOTCANDS>){
	chomp;
	$pulsefile = $_;
	@pulsedata = split(/\_/,$pulsefile);
	$snr = $pulsedata[1];
	$start = $pulsedata[2];
	$width = $pulsedata[4];
	$dec = $pulsedata[5];
	$unformatteddm = $pulsedata[6];
	$dm = sprintf "%07.2f",$pulsedata[6];
	$maxbeam = $pulsedata[-1];
	$maxbeam =~ s/.PULSE//;
	$maxbeam =~ s/.pulse//;
	$plotnamebase = $basename."_at_".$start."_DM_".$dm."_SNR_".$snr."beam$maxbeam";
    

	$ib = $maxbeam;
#	foreach $ib (@ibeam){
	$ibchar = convertib($ib);
	$datafilebase = $origbasename;
	system "ln -s $datafilebase.$dm.tim $plotnamebase.beam$ibchar";
	# Set up links for beam of britest detection
	if ($ib == $maxbeam || $ibchar eq $maxbeam){
	    print STDERR "Maxbeam is $maxbeam\nTrolling file $datafilebase\n\n";
	    system "ln -s $infile $plotnamebase.fil";
	    system "ln -s $datafilebase.0000.00.tim $plotnamebase.dm0";
		
	    # Calculate delays for plotting later
	    $sampdelay = `dmdelay $infile $dm`;
	    $pulsespan = $sampdelay + $width;
	    if ($start-$pulsespan < 0){ $sampskip = 0;}
	    else { $sampskip = $start - $pulsespan; }
	    $sampread = 3*$pulsespan;
	    $sampread = int($sampread+0.5);
	    $sampskip = int($sampskip+0.5);

	    # If candidate has a sweep smaller than one scrunched
	    # sample or there are only a few points on the plot,
	    # lower the decimation to allow better viewing
	    if ($sampread/$dec<6 && $dec>1){
		$dec = $dec/2;
		print STDERR "\n\tWARNING:\n\tReducing decimation factor to make plot look good.\n\n";
		system "sleep 0.5";
	    }
	    
	    $i = 0;
	    foreach $arrdm (@dmlist){
		if ($arrdm == $dm){ last; }
		else { $i++; }
	    }
	    system "ln -s $pulsefile $plotnamebase.PULSEDATA";
	    system "ln -s $datafilebase.$dmlist[$i-1].tim $plotnamebase.dmlo";
	    if ($dmlist[$i+1] ne ''){system "ln -s $datafilebase.$dmlist[$i+1].tim $plotnamebase.dmhi";}
	} else {
	    print STDERR "\n\tERROR:\n\tThere is something wrong with the beam names.\n\n";
	    exit(6);
	}

	$beampad = $maxbeam;
	if ($beampad eq 'A'){
	    $beampad = 10;
	} elsif ($beampad eq 'B'){
	    $beampad = 11;
	} elsif ($beampad eq 'C'){
	    $beampad = 12;
	} elsif ($beampad eq 'D'){
	    $beampad = 13;
	}
	
	$beampad = '0'.$beampad;
	$beampad = substr($beampad,-2);
	print STDERR "------------------------------------\nPlotting $plotnamebase:\n";
	print STDERR "\n\nquickgplot $plotnamebase -s $sampskip -r $sampread -dec $dec -bestbeam $maxbeam";
	if ($killmask eq ''){
	    system "quickgplot $plotnamebase -s $sampskip -r $sampread -dec $dec -bestbeam $maxbeam\n\n";
	} else {
	    system "quickgplot $plotnamebase -s $sampskip -r $sampread -dec $dec -bestbeam $maxbeam -k $killmask\n\n";
	}
	system "rm $plotnamebase.dm0 $plotnamebase.dm?? $plotnamebase.beam? $plotnamebase.PULSEDATA $plotnamebase.fil";
    }
    $nplots = `ls -1 $basename*.png | wc -l`;
    chomp($nplots);
    if ($nplots == 0){
	print STDERR "\n\tERROR:\n\tSomething might be wrong with your plotting software...\n\n";
    } elsif ($nplots != $ncands){
	print STDERR "\n\tWARNING:\n\tWas only able to make plots for $nplots out of $ncands candidates.\n\n";
    } else {
	print STDERR "\n\n\t$nplots candidates plotted successfully!\n\n";
    }
    if (!-e "tim-files-$basename"){
	system "mkdir tim-files-$basename";
    }
    system "mv $basename*.tim tim-files-$basename";
}




sub convertib(){
    my $beam = shift @_;
    my $charbeam = 0;
    if ($beam<10){
	$charbeam = int($beam);
    } elsif ($beam==10){
	$charbeam = "A";
    } elsif ($beam==11){
	$charbeam = "B";
    } elsif ($beam==12){
	$charbeam = "C";
    }elsif ($beam==13){
	$charbeam = "D";
    }
    return($charbeam);
}


sub findindex(){
    my ($testval,@arr) = @_;
    my $index = 0;
    my $i;
    for ($i=0;$i<@arr;$i++){
	if ($testval eq $arr[$i]){
	    $index = $i;
	    return($index);
	}
    }
    return(-1);
}



sub printhelp{
print STDERR "

***************************************
Large Burst Searching Program Help Menu
***************************************

USAGE:
\tspsearch <filename(s)> [-<options>]

To use this program correctly, one MUST specify at least one file name, and
optionally other parameters as given below (default values when parameter not
set by user are given in brackets).\n"; #The script will operate on a copy of the input file that is called, as long as the working directory is not set as the same location as the data."

print STDERR "
\t-box\t\tBoxcar-search up to a filter width of N samples (must be power of 2) [512]
\t-k\t\tInput channel mask [no channel mask]
\t-mindm\t\tLowest DM to search (program will calculate steps) [0]
\t-maxdm\t\tHighest DM to search (program will calculate steps) [500]
\t-sig\t\tSearch for pulses above this significance [6]
\t-smear\t\tDM smearing tolerance (influences # DM trials; e.g. lower tolerance, more trials) [1.25]
\t-wid\t\tPulse separation tolerance (units: # bins) [30]
\t-zapdm\t\tDM below which considered RFI if pulse peaks there in pc/cc [default: 3]

"
}
