#include <cstdlib>
#include <vector>
#include <cmath>
#include <stdint.h>
#include "mjklog.h"
#include "mjk_cmd.h"

#include "sigproc/filfile.hpp"
#include "sigproc/outputfilfile.hpp"
#include "sigproc/filterbankblock.hpp"

void merge_blocks(std::vector<sigproc::FilterbankBlock*> input_blocks, sigproc::FilterbankBlock *output_block);

void print_usage(){
    printf("Merge .fil files\n");
    printf("================\n");
    printf("\n");
    printf("Take many fil files, sensibly combine them. They must align in time and be consecutive or align in frequency.\n");
    printf("--help,-h           : print this help text\n");
    printf("\n");
}

int main(int argc, char** argv) {


    if(argc<2 || getB("--help","-h", argc,argv, false)){
        print_usage();
        exit(0);
    }


    char *outfile = getS("--output","-o", argc,argv,"out.fil");

    std::vector<sigproc::FilFile> input_files;

    getArgs(&argc, argv);
    if (argc < 2){
        logerr("No files specified on command line...");
        print_usage();
        exit(1);
    }

    for(int iarg=1; iarg < argc; ++iarg) {
        sigproc::FilFile infile(argv[iarg]);
        infile.initialise();
        infile.debug();
        if (infile.valid()) {
            logmsg("loaded '%s' ok.",argv[iarg]);
            input_files.push_back(infile);
        } else {
            logerr("Failed to load '%s'.",argv[iarg]);
        }
    }

    sigproc::OutputFilFile output(outfile,input_files[0]);
    for (int i=0; i < input_files.size(); ++i){
        bool good = output.expandTo(input_files[i]);
        if (!good){
            logerr("Cannot merge %s",input_files[i].filename());
        }
    }
    output.setNbits(32);
    output.debug();
    logmsg("Initialise");
    output.initialise();

    int blocklength = 8192;

    uint64_t sample = 0;
    while (!input_files[0].eof()){
        logmsg("%lf s",sample*output.sample_interval());
        std::vector<sigproc::FilterbankBlock*> inputblocks;
        for (int i=0; i < input_files.size(); ++i){
            logmsg("Readblock: %s",input_files[i].filename());
            sigproc::FilterbankBlock* block = input_files[i].readBlock(sample,blocklength);
//            logmsg("Push block");
            inputblocks.push_back(block);
//            logmsg("Pushed block");
        }
//        logmsg("make output block");
        sigproc::FilterbankBlock* outblock = output.createBlock(blocklength);
//        logmsg("merge blocks");
        merge_blocks(inputblocks,outblock);
//        logmsg("write blocks");
        output.writeBlock(outblock);
        //logmsg("zz %lg %lg",inputblocks[0]->_data[0],outblock->_data[0]);
        sample += blocklength;
        delete outblock;
        for (int i=0; i < input_files.size(); ++i){
            delete (inputblocks[i]);
        }
    }


    return 0;
}

void merge_blocks(std::vector<sigproc::FilterbankBlock*> input_blocks, sigproc::FilterbankBlock *output_block){
    float* weights=static_cast<float*>(calloc(sizeof(float),output_block->_raw_length));
    int nifs=output_block->_nifs;
    for (int ib=0; ib < input_blocks.size(); ++ib){
        int chanoffset = output_block->_filfile->getChanOffset(*input_blocks[ib]->_filfile);
//        logmsg("Chan offset %d",chanoffset);
        for(int isamp=0; isamp < input_blocks[ib]->_length; ++isamp){
            float* out = output_block->_data + output_block->_nchans*nifs*isamp + chanoffset;
            float* in = input_blocks[ib]->_data + input_blocks[ib]->_nchans*nifs*isamp;
            //logmsg(" isamp: %d  %d   %d   %d %d",isamp,in-input_blocks[ib]->_data,out-output_block->_data,nifs,input_blocks[ib]->_nchans);
            for(int ich=0; ich < input_blocks[ib]->_nchans; ++ich){
                out[ich] += in[ich];
            }
        }
    }

}

