#pragma once
#include <cstdio>
#include <stdint.h>


namespace sigproc {
class OutputFilFile;
class FilterbankBlock;

class FilFile {

    public:

        /**
         * Default constructor, takes a filename
         */
        FilFile(const char* filename);

        virtual void initialise();

        FilterbankBlock *readBlock(int start_sample, int length);

        bool eof() const;

        int getChanOffset(const FilFile &other) const;
        /**
         * Is the file readable?
         */
        bool valid() const;

        void debug() const;

        const char* filename() const{
            return _filename;
        }

        double sample_interval() const {
            return _tsamp;
        }

    protected:

        /**
         * Get the sigproc global parameters into this file
         */
        void globalToThis();

        void thisToGlobal();

        char _source_name[80];
        int _machine_id, _telescope_id, _data_type, _nchans, _nbits, _nifs, _scan_number,
              _barycentric,_pulsarcentric; /* these two added Aug 20, 2004 DRL */
        double _tstart,_mjdobs,_tsamp,_fch1,_foff,_refdm,_az_start,_za_start,_src_raj,_src_dej;
        char _isign;

        bool _valid;
        int _header_length;

        const char* _filename;
        FILE* _rawfile;


        friend class OutputFilFile;
        friend class FilterbankBlock;

} ; // class Filfile
} // namespace sigproc

