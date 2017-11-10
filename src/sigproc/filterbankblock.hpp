#pragma once
#include <sigproc/filfile.hpp>
#include <stdint.h>


namespace sigproc {
class FilterbankBlock {
    public:
        FilterbankBlock(uint64_t start, uint64_t length, const FilFile *filfile);
        ~FilterbankBlock();
        float* _data;

        const int _nchans;
        const int _nifs;

        const uint64_t _start;
        const uint64_t _length;
        const uint64_t _raw_length;
        const FilFile *_filfile;
    private:
}; //class FilterbankBlock

} // namespace sigproc


