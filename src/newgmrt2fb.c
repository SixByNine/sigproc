#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#include "filterbank.h"
void newgmrt2fb(FILE *input, FILE *output) /* includefile*/
{
    unsigned char *buffer = (unsigned char*)malloc(nchans);
    int ich;
    int flipch;
    unsigned char temp;

    if (nbits==8){
        while (!feof(input)){
            fread(buffer,1,nchans,input);
            for(ich=0; ich< nchans; ich++){
                buffer[ich] &= gmrtzap[ich];
            }
            if(invert_band){
                for(ich=0; ich < nchans/2; ich++){
                    flipch = nchans - ich-1;
                    SWAP(buffer[ich],buffer[flipch]);
                }
            }

            fwrite(buffer,1,nchans,output);

        }
    } else if(nbits==4) {
        unsigned char* packzap = (unsigned char*)calloc(nchans,1);
        for(ich=0; ich < nchans; ich++) packzap[ich]=0xFF;

        for(ich=0; ich < nchans; ich++){
            if(ich%2){
                // even chans
                packzap[ich] &= gmrtzap[ich]>>4;
            } else {
                packzap[ich] &= (gmrtzap[ich]&0x0F)<<4;
            }
        }

        while (!feof(input)){
            fread(buffer,1,nchans/2,input);
            for(ich=0; ich< nchans; ich++){
                buffer[ich/2] &= packzap[ich];
            }

            if (invert_band && nbits==4){
                // re-order bytes
                for(ich=0; ich < nchans/4; ich++){
                    flipch = nchans/2 - ich-1;
                    SWAP(buffer[ich],buffer[flipch]);
                }
                // also flip all the pairs of channels...
                // aaaabbbb => bbbbaaaa
                for(ich=0; ich < nchans/2; ich++){
                    char a = buffer[ich]>>4;
                    char b = (buffer[ich]&0x0F) << 4;
                    buffer[ich] = a&b;
                }
            }
            fwrite(buffer,1,nchans/2,output);
        }

    }

}
