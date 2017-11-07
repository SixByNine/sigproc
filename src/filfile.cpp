
#include <cstring>

#include "sigproc/filfile.hpp"

#include "header.h"
#include "sigproc.h"
#include "mjklog.h"



sigproc::FilFile::FilFile(const char* filename) : _valid(false), _filename(filename) {
}

void sigproc::FilFile::initialise() {
    _rawfile = fopen(_filename,"rb");
    if (_rawfile){
        _header_length = read_header(_rawfile);
        if (data_type == 0 || data_type == 1 ){
            globalToThis();
            _valid = true;
        }
    }
}


bool sigproc::FilFile::valid() const {
    return _valid;
}

void sigproc::FilFile::debug() const {
    logmsg("'%s' '%s' %lf %lf %ld",_filename,_source_name,_fch1,_foff,_nchans);
}


void sigproc::FilFile::globalToThis() {
    strncpy(_source_name,source_name,80);
    _machine_id   = machine_id;
    _telescope_id = telescope_id;
    _data_type    = data_type;
    _nchans       = nchans;
    _nbits        = nbits;
    _nifs         = nifs;
    _scan_number  = scan_number;
    _barycentric  = barycentric;
    _pulsarcentric = pulsarcentric;
    _tstart       = tstart;
    _mjdobs       = mjdobs;
    _tsamp        = tsamp;
    _fch1         = fch1;
    _foff         = foff;
    _refdm        = refdm;
    _az_start     = az_start;
    _za_start     = za_start;
    _src_raj      = src_raj;
    _src_dej      = src_dej;
    _nbits        = nbits; 
    _isign        = isign;
}


