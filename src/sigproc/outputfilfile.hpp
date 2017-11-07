#pragma once
#include <sigproc/filfile.hpp>

namespace sigproc {
class OutputFilFile : public FilFile {
    public:
        OutputFilFile(const char* filename);
        OutputFilFile(const char* filename,const FilFile& copy);

        virtual void initialise();

        bool expandTo(const FilFile& other);

    private:

}; // class OutputFilFile
}
