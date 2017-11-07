#pragma once
#include <cstdio>


namespace sigproc {
class OutputFilFile;

class FilFile {

    public:

        /**
         * Default constructor, takes a filename
         */
        FilFile(const char* filename);

        virtual void initialise();


        /**
         * Is the file readable?
         */
        bool valid() const;

        void debug() const;

        const char* filename() const{
            return _filename;
        }

    protected:

        /**
         * Get the sigproc global parameters into this file
         */
        void globalToThis();

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

} ; // class Filfile
} // namespace sigproc

