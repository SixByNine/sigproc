
#include <sigproc/outputfilfile.hpp>
#include <cmath>

#include "mjklog.h"



sigproc::OutputFilFile::OutputFilFile(const char* filename) : FilFile(filename) {
}


sigproc::OutputFilFile::OutputFilFile(const char* filename, const FilFile &copy) : FilFile(copy) {
    _filename = filename;
}


void sigproc::OutputFilFile::initialise() {
}

bool sigproc::OutputFilFile::expandTo(const FilFile &other) {
    if (other._foff != _foff) return false;
    if (other._tsamp != _tsamp) return false;
    if (other._foff * _foff < 0) return false;

    double start_chan = (other._fch1 - _fch1) / _foff;
    if (fabs(start_chan - round(start_chan)) > 0.01){
        return false;
    }
    
    int istart_chan = (int)(round(start_chan));

    if (start_chan < 0){
        // we need to start at an earlier channel.
        _fch1 += start_chan * _foff;
        _nchans -= start_chan;
    }


    double end_chan = (other._fch1+other._nchans*other._foff - _fch1) / _foff;
    if (fabs(end_chan - round(end_chan)) > 0.01){
        return false;
    }
    
    int iend_chan = (int)(round(end_chan));
    if (iend_chan > _nchans) {
        _nchans = iend_chan;
    }

    return true;
}
