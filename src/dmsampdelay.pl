#!/usr/bin/perl
#
#  Dispersion delay calculator.
#  Sarah B. Spolaor, Swinburne University
#  (c) 2007
#
#  Given an input filterbank file and dispersion measure (units pc/cc), this
#  script returns the number of filterbank samples for a pulse of the given DM
#  to cross the file's bandwidth.
#

$file = $ARGV[0];
$dm = $ARGV[1];

if ($#ARGV != 1){
    print STDERR "\n\tINPUT ERROR in dmdelay:\n\tINPUTS: <file (usually a .fil file)> <dm (pc/cc)>\n\n";
    exit(1);
} elsif (!-e $file){
    print STDERR "\n\tERROR in dmdelay:\n\tInput file $file does not exist.\n\n";
    exit(2);
} else{
    open(FILE, "header $file |");
    while(<FILE>){
	chomp;
	@array = split(/\s+/,$_);
	$nametest = $array[0].$array[1];
	if ($nametest eq 'Sampletime'){$tsamp = $array[4]/1000;} #sample time in milliseconds
	if ($nametest eq 'Frequencyof'){$f1 = $array[6];}
	if ($nametest eq 'Channelbandwidth'){$bw = $array[4]};
	if ($nametest.$array[2] eq 'Numberofchannels') {$nchan = $array[4]};
    }
    close(FILE);
    system "rm -f $tempfile";
    $f2 = $f1 + $nchan * $bw;
    $delayinms = $dm * 4.15 * (($f2/1000)**(-2)-($f1/1000)**(-2));
    $delayinsamples = $delayinms / ($tsamp);
    printf "%d",$delayinsamples+1;
}
