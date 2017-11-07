#include <cstdlib>
#include <vector>
#include "mjklog.h"
#include "mjk_cmd.h"

#include "sigproc/filfile.hpp"
#include "sigproc/outputfilfile.hpp"

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
    output.debug();
    output.initialise();



    return 0;
}
