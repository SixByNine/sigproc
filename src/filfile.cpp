
#include <cstring>
#include <cmath>

#include "sigproc/filfile.hpp"
#include "sigproc/filterbankblock.hpp"

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

bool sigproc::FilFile::eof() const {
    return feof(_rawfile);
}

int sigproc::FilFile::getChanOffset(const sigproc::FilFile &other) const {
    double start_chan = (other._fch1 - _fch1) / _foff;
    int istart_chan = (int)(round(start_chan));
    return istart_chan;
}



sigproc::FilterbankBlock *sigproc::FilFile::readBlock(int start_sample, int length) {
    FilterbankBlock *block = new FilterbankBlock(start_sample,length,this);

    // get to the right place in the file.
    //fseek(_rawfile,_header_length+(start_sample*_nbits*_nchans*_nifs)/8,SEEK_SET);
//int read_block(FILE *input, int nbits, float *block, int nread) /*includefile*/
    read_block(_rawfile, _nbits,block->_data, length*_nchans*_nifs);

    return block;
}


void sigproc::FilFile::debug() const {
    logmsg("'%s' '%s' %lf %lf %ld %d nbits:%d",_filename,_source_name,_fch1,_foff,_nchans,_nifs,_nbits);
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

void sigproc::FilFile::thisToGlobal() {
    strncpy(source_name,_source_name,80);
    machine_id   = _machine_id;
    telescope_id = _telescope_id;
    data_type    = _data_type;
    nchans       = _nchans;
    nbits        = _nbits;
    nifs         = _nifs;
    scan_number  = _scan_number;
    barycentric  = _barycentric;
    pulsarcentric = _pulsarcentric;
    tstart       = _tstart;
    mjdobs       = _mjdobs;
    tsamp        = _tsamp;
    fch1         = _fch1;
    foff         = _foff;
    refdm        = _refdm;
    az_start     = _az_start;
    za_start     = _za_start;
    src_raj      = _src_raj;
    src_dej      = _src_dej;
    nbits        = _nbits;
    isign        = _isign;
}

