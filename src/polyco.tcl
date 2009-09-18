###############################################################################
# POLYCO: A Tcl script to generate a polyco.dat file by running TEMPO
# Modification history...
# 02/01/30 - changed to tcl script rather than expect (not always available)
# 05/06/24 - added par option and strip of J in front of psrname if need be
###############################################################################
set timeout -1

if { ! [ info exists env(TEMPO) ] } {
  set env(TEMPO) /home/pulsar/tempo11
}

proc help {} {
    global freq nspan ncoeff maxha site
    puts ""
    puts "polyco - a script to run TEMPO to generate a polyco.dat file"
    puts ""
    puts "usage: polyco psrname -{options}"
    puts ""
    puts "psrname   - name of the pulsar as it appears in tztot.dat"
    puts ""
    puts "-freq   f - frequency in MHz (def=$freq MHz)"
    puts "-nspan  n - span of each polyco set in minutes (def=$nspan min)"
    puts "-ncoeff n - number of coefficients in each polyco set (def=$ncoeff)"
    puts "-maxha  h - maximum hour angle (def=$maxha hours)"
    puts "-mjd    m - mjd to calculate for (def=today)"
    puts "-mjds   s - starting mjd (def=today)"
    puts "-mjdf   f - finishing mjd (def=today)"
    puts "-par    p - specify parfile (def=tzpar area)"
    puts "-site   s - specify site code or alias (def=$site)"
    puts ""
}
set value 0
set freq 1410
set nspan 15
set ncoeff 9
set maxha 2
set site Arecibo
set par ""
if {$argv == "" || $argv == "help"} {
    help
    exit
}
set psrname [string trimleft [lindex $argv 0] J]
set mjds [set mjdf ""]
foreach item [lrange $argv 1 end] {
    if {!$value} {
	set item [string range $item 1 end]
	switch $item {
	    freq     -
	    nspan    -
	    ncoeff   -
	    maxha    -
	    mjd      -
	    mjds     -
	    site     -
	    par      -
	    mjdf     {set value 1; set key $item}
	    default  {help;puts "invalid key $item";exit}
	}	
    } else {
	switch $key {
	    freq      {set freq $item}
	    nspan     {set nspan $item}
	    ncoeff    {set ncoeff $item}
	    maxha     {set maxha $item}
	    mjds      {set mjds $item}
	    mjdf      {set mjdf $item}
	    site      {set site $item}
	    par       {if [file exists $item] {set par  "-f $item"} else {puts "parameter file $item not found";exit}}
	    mjd       {set mjds [expr $item-1]; set mjdf [expr $item+1]}
	}
	set value 0
    }
}
if {$mjds != "" && $mjdf == ""} {set mjdf $mjds}
if {$mjdf != "" && $mjds == ""} {set mjds $mjdf}
switch [string tolower $site] {
    gb -
    gbt -
    1 {set site 1}
    quabbin -
    2 {set site 2}
    ao - 
    arecibo -
    3 {set site 3}
    hobart -
    4 {set site 4}
    princeton -
    5 {set site 5}
    vla -
    6 {set site 6}
    parkes -
    7 {set site 7}
    lovell -
    jodrell -
    jb -
    8 {set site 8}
    torun -
    torun32m -
    9 {set site 9}
    gb140ft -
    a {set site a}
    gb85ft -
    b {set site b}
    vlasite -
    c {set site c}
    bologna -
    d {set site d}
    most -
    e {set site e}
    nancay -
    f {set site f}
    effelsberg -
    g {set site g}
    kalyazin -
    h {set site h}
    fallbrook - 
    i {set site i}
    gmrt - 
    r {set site r}
    default {
	help
	puts "ERROR: site $site unknown!"
	set obsysfile $env(TEMPO)/obsys.dat
	if [file exists $obsysfile] {
	    puts "consult $obsysfile"
	}
	exit
    }
}
set file [open polyco.dat w]
close $file
set file [open tz.in w]
puts $file "$site    2   30   9  1410" 
puts $file ""
puts $file ""
puts $file "$psrname $nspan $ncoeff $maxha $freq"
close $file
set file [open tempo.old w]
puts $file  "#!/bin/csh"
puts $file  "tempo -z $par << END" 
puts $file  "$mjds $mjdf"      
puts $file  "END"      
close $file
exec chmod +x tempo.old
catch {exec ./tempo.old}
if ![file size polyco.dat] {
    set file [open tempo.new w]
    puts $file  "#!/bin/csh"
    puts $file  "tempo -z << END" 
    puts $file  "$mjds $mjdf"      
    puts $file  "END"      
    close $file
    exec chmod +x tempo.new
    catch {exec ./tempo.new}
    if ![file size polyco.dat] {
	puts "Error running TEMPO!"
	puts "check files tz.in tempo.old tempo.new..."
	exit
    } else {
	exec rm -f tempo.new
    }
}
puts "TEMPO appears to have run successfully..."
puts [exec ls -l polyco.dat]
exec rm -f tz.in tempo.old
exit
