
#include <assert.h>
#include <string>
#include <cstdio>
#include <cmath>
#include <iostream>


#include <sigproc.h>
#include <header.h>

extern "C" {
    char ifstream[8];
    int obits;
}

void rescale_and_trim(float* inbuf, unsigned char* outbuf, double* scale, double* offset, int raw_nchan, int out_nchan, int skipchan, int nsamp);
int read_fil_to_buf(char *filename, float *&buf);


/**
 * PAF redigitisation, rechannelisation and merging code
 * M. Keith & M. Malenta 2018
 */
int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Error - need at least two arguments" << std::endl <<std::endl;
        std::cout << "PAF redigitisation, rechannelisation and merging code" << std::endl;
        std::cout << "M. Keith & M. Malenta 2018" << std::endl;
        std::cout << "This program is only useful for working on the paf data from" << std::endl;
        std::cout << "M.Malenta's filterbank code." << std::endl << std::endl;
        std::cout << "Example Usage:" << std::endl;
        std::cout << argv[0] << " out.fil in*.fil" << std::endl;
        exit(1);
    }


    assert(argc > 2);

    std::cout << "Output to: " << argv[1] << std::endl;
    std::cout << "Reading:" << std::endl;
    std::cout << argv[2] << std::endl;


    // these parameters describe the input and output data.
    const unsigned int raw_nchan=567;
    const unsigned int output_nchan=512;
    const unsigned int skip_chans=27;

    float *inbuf;
    // pointer is passed by reference here. We have to delete it later.
    int nsamp = read_fil_to_buf(argv[2],inbuf);


    // these arrays are used to store the per-channel mean and variance.
    double means[output_nchan];
    double var[output_nchan];


    // I am not sure that this is nessicary in C++, but better to be safe.
    for (int ichan = 0; ichan < output_nchan ; ++ichan) {
        means[ichan]=0;
        var[ichan]=0;
    }


    // I know that a two-pass mean/variance is slow, but it is reliable and we have
    // cpu cycles to spare. We only compute this once per observation chunk
    for (int isamp = 0; isamp < nsamp ; ++isamp) {
        // this points to the start of the channels we want to output
        float *sbuf = inbuf + isamp*raw_nchan + skip_chans;
        for (int ichan = 0; ichan < output_nchan ; ++ichan) {
            means[ichan] += sbuf[ichan];
        }
    }

    for (int ichan = 0; ichan < output_nchan ; ++ichan) {
        means[ichan] /= static_cast<double>(nsamp);
    }

    for (int isamp = 0; isamp < nsamp ; ++isamp) {
        // again, this points to the start of the channels we want to output
        float *sbuf = inbuf + isamp*raw_nchan + skip_chans;
        for (int ichan = 0; ichan < output_nchan ; ++ichan) {
            double val = sbuf[ichan]-means[ichan];
            var[ichan]   += val*val;
        }
    }

    for (int ichan = 0; ichan < output_nchan ; ++ichan) {
        var[ichan] /= static_cast<double>(nsamp);
    }




    // open the output file, and also a file to store the weights used.
    char* outfilname=argv[1];
    FILE* outfilptr = fopen(outfilname,"w");
    std::string logfname(outfilname);
    logfname += ".scale"; // weights file is outfile.fil.scale
    FILE* logfile = fopen(logfname.c_str(), "w");

    double offset[nchans];
    double scale[nchans];
    for (int ichan = 0; ichan < output_nchan ; ++ichan) {
        offset[ichan] = means[ichan];
        if (var[ichan] == 0) {
            scale[ichan] = 1.0;
        } else {
            scale[ichan]  = 24.0/(sqrt(var[ichan]));
        }
        if (std::isnan(scale[ichan])) {
            std::cout << "NAN "<< ichan << std::endl;
        }
        fprintf(logfile,"%03d % 14.8lg % 14.8lg % 14.8lg\n",ichan,offset[ichan],scale[ichan],var[ichan]);
    }
    fclose(logfile);

    // this will store the output buffer
    unsigned char *outbuf = new unsigned char[output_nchan*nsamp];
    rescale_and_trim(inbuf, outbuf, scale, offset, raw_nchan, output_nchan, skip_chans,nsamp);


    // adjust the global variables used by sigproc
    fch1 += skip_chans*foff;
    nchans = output_nchan;
    obits=8;
    ifstream[0]='Y';
    machine_id=88; // for now hard code this to the PAF number

    // this writes the header
    filterbank_header(outfilptr);
    fwrite(outbuf,1,nchans*nsamp,outfilptr);

    delete[] inbuf;
    delete[] outbuf;
    // 0 is prog name
    // 1 is output file
    // 2 is file we already did
    for (int ifile = 3; ifile < argc; ifile++){
        std::cout << argv[ifile] << std::endl;
        // read the rest of the data...
        int nsamp = read_fil_to_buf(argv[ifile],inbuf);
        unsigned char *outbuf = new unsigned char[output_nchan*nsamp];
        rescale_and_trim(inbuf, outbuf, scale, offset, raw_nchan, output_nchan, skip_chans,nsamp);
        fwrite(outbuf,1,output_nchan*nsamp,outfilptr);
        delete[] inbuf;
        delete[] outbuf;
    }

    fclose(outfilptr);
    return 0;
}

void rescale_and_trim(float* inbuf, unsigned char* outbuf, double* scale, double* offset, int raw_nchan, int out_nchan, int skipchan, int nsamp) {
    for (int isamp = 0; isamp < nsamp ; ++isamp) {
        float *sbuf = inbuf + isamp*raw_nchan + skipchan;
        unsigned char *osbuf = outbuf + isamp*out_nchan;
        for (int ichan = 0; ichan < out_nchan ; ++ichan) {
            double v = (sbuf[ichan]-offset[ichan])*scale[ichan] + 64;
            if (v > 255.) v=255.;
            if (v < 0.) v=0.;
            osbuf[ichan] = (unsigned char)v;
        }
    }
}

int read_fil_to_buf(char *filename, float *&buf){
    FILE *file_ptr = fopen(filename, "r");

    int hdrlen = read_header(file_ptr);
    int nsamp = nsamples(filename,hdrlen, nbits,nifs,nchans);

    assert(nbits==32);
    assert(nifs==1);
    assert(nchans==567);

    buf = new float[nsamp*nchans];
    int nread = fread(buf, 4, nsamp*nchans, file_ptr);

    assert(nread==nsamp*nchans);

    return nsamp;
}


