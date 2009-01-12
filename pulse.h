typedef struct plinfo {
  int index;                          /* the position of the pulse         */
  float amp;                          /* the amplitude of the pulse in     */
                                      /* of the standard deviation         */
  double mean;                        /* the mean of the "chunk" where the */
                                      /* pulse was found                   */
  double rms;                         /* the rms of the "chunk" where the  */

} Pulsus;

int iindx;                            /* Since we don't read in the entire */
                                      /* file at once, we must keep track  */
                                      /* of how many reads have been done  */
                                      /* on the file in order to get the   */
                                      /* pulse index correct in the time   */
                                      /* series.                           */

char *optarg;
int opterr;
int optind;

unsigned short iterate;               /* Use iterative method for mean     */
                                      /* and rms?  iterate==0 means no     */
float thresh;                         /* The threshold in units of the rms */

char file_name[80];                   /* Basename of input file            */
char input_file[80];                  /* Full name of input file           */

int ndm;                              /* dm channel number from file name  */ 

int nsubband;                         /* subband number from file name     */

int nsmax;                            /* maximum number of smoothings to try */















