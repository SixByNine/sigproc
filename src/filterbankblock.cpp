#include "sigproc/filterbankblock.hpp"
#include "mjklog.h"
#include <cstdlib>


sigproc::FilterbankBlock::FilterbankBlock(uint64_t start, uint64_t length, const FilFile *filfile) :
    _start(start), _length(length), _filfile(filfile), _nchans(filfile->_nchans),_nifs(filfile->_nifs),_raw_length(length*_nchans*_nifs) {
        _data = static_cast<float*>(calloc(sizeof(float),_raw_length));
      //  logmsg("Allocate block %lp",_data);
    }

sigproc::FilterbankBlock::~FilterbankBlock() {
    //logmsg("deallocate block %lp",_data);
    free(_data);
}

