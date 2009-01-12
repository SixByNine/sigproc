#!/bin/csh
# This script updates the include file sigproc.h for the sigproc library
set v = `awk '{print $1}' version.history | tail -2 | head -1`
cat << END1   >! sigproc.h
/* sigproc.h: Automatically generated include file for sigproc-$v */
#ifdef __cplusplus
extern "C" {
#endif
#include "polyco.h"
#include "epn.h"
#include "version.h"
END1
grep includefile *.c | awk -F"/*" '{print $1";"}' \
	| awk -F: '{print $2}' | sort >> sigproc.h
cat << END2   >> sigproc.h
#ifdef __cplusplus
}
#endif
END2
chmod 666 sigproc.h
exit
