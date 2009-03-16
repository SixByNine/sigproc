###############################################################################
# MAKEDMLIST - a tcl script to make a list of DMs given a filterbank file.....
###############################################################################
set file [lindex $argv 0]
if {[llength $argv] > 1} {
	set dmmin [lindex $argv 1]
} else {
	set dmmin 0.0
}
if {[llength $argv] > 2} {
	set dmmax [lindex $argv 2]
} else {
	set dmmax 100.0
}
if ![file exists $file] {
	puts "usage: makedmlist name_of_filterbank_file (dmmin) (dmmax)"
	puts "DM minimum and maximum ranges are optional (0-100 is default)"
	exit
}
set tsamp [expr [exec header $file -tsamp]/1000.0]
set nchan [exec header $file -nchans]
set fcent [exec header $file -fmid]
set bw    [exec header $file -bandwidth]
set j 0
set dm 0.0
while {1} {
  for {set i 1} {[expr $i<=$nchan]} {incr i} {
        set n [expr $i.0-1.0]
        set dm [format %.1f [expr 1.205e-7*$n*$tsamp*pow($fcent,3.0)/$bw]]
        if {[expr $dm > $dmmax]} {
                foreach dm [lsort -real -increasing -unique $dmlist] {
                        if {$dm > $dmmin} {puts $dm}
                }
                exit
        }
        lappend dmlist $dm
  }
  set tsamp [expr $tsamp*2.0]
  incr j
}
