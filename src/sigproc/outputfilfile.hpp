#pragma once
#include <sigproc/filfile.hpp>

namespace sigproc {
class OutputFilFile : public FilFile {
    public:
        OutputFilFile(const char* filename);
        OutputFilFile(const char* filename,const FilFile& copy);

        virtual void initialise();

        void writeBlock(const FilterbankBlock *block);
        FilterbankBlock *createBlock(int length);

        void setNbits(int nb) {
            _nbits=nb;
        }

        bool expandTo(const FilFile& other);

    private:

        int _current_sample;

}; // class OutputFilFile
}
