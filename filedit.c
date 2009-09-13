#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

int read_header(FILE *inputfile);

void fix_header(FILE* file, char* newname, double newra, double newdec);
void zap_em(FILE* file, int* tzaps[2], int ntzaps, int fzaps[1024], int nfzaps,float mean, float sigma);
float get_random_value(float mean, float sigma);

void print_usage(){
	printf("Modify a .fil file\n");
	printf("==================\n");
	printf("\n");
	printf("Modifies a .fil file 'in place'\n");
	printf("Options:\n");
	printf("--ra,-r {ra}         : modify the ra to {ra}. in form hhmmss.xxx \n");
	printf("--dec,-d {dec}       : modify the dec to {dec}. in form ddmmss.xxx \n");
	printf("--src-name,-n        : modify the source name.\n");
	printf("--time-zap,-t \"s e\"  : 'zap' samples between 's' and 'e'\n");
	printf("\n");

}

int main (int argc, char** argv){

	struct option long_opt[256];
	const char* args = "d:hn:t:r:m:s:k:";
	int opt_flag = 0;
	int long_opt_idx = 0;
	int* timezaps[2];
	int fzaps[1024];
	int ntimezaps,nfreqzaps;
	double newra;
	double newdec;
	float mean,sigma;
	char newname[1024];
	int c;
	FILE* file;
	char killfile[20];
	int killswitch=0;
	int arr_size = 1024;

	newname[0]='\0';
	ntimezaps=0;
	nfreqzaps=0;
	mean=0;
	sigma=1;
	
	newra =900000000;
	newdec=900000000;


	long_opt[long_opt_idx].name = "help";
	long_opt[long_opt_idx].has_arg = no_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'h';

	long_opt[long_opt_idx].name = "time-zap";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 't';

	long_opt[long_opt_idx].name = "mean";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'm';

	long_opt[long_opt_idx].name = "sigma";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 's';


	long_opt[long_opt_idx].name = "ra";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'r';

	long_opt[long_opt_idx].name = "dec";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'd';

	long_opt[long_opt_idx].name = "src-name";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'n';

	long_opt[long_opt_idx].name = "tkill";
	long_opt[long_opt_idx].has_arg = required_argument;
	long_opt[long_opt_idx].flag = NULL;
	long_opt[long_opt_idx++].val = 'k';


	long_opt[long_opt_idx].name = 0;
	long_opt[long_opt_idx].has_arg = 0;
	long_opt[long_opt_idx].flag = 0;
	long_opt[long_opt_idx++].val = 0;


	while ((c = getopt_long(argc, argv, args, long_opt, &long_opt_idx)) != -1) {
		switch (opt_flag) {
			default:
				opt_flag = 0;
				switch (c) {
					case 'h':
						print_usage();
						exit(0);
						break;
					case 't':
						sscanf(optarg,"%d %d",timezaps[0]+ntimezaps,timezaps[1]+ntimezaps);
						ntimezaps++;
						break;
					case 'r':
						newra = atof(optarg);
						break;
					case 'd':
						newdec=atof(optarg);
						break;
					case 'n':
						strcpy(newname,optarg);
						break;
					case 'm':
						mean=atof(optarg);
						break;
					case 's':
						sigma=atof(optarg);
						break;
				        case 'k':
					  sscanf(optarg, "%s",killfile);
					  killswitch = 1;
					        break;
				}
		}
	}

	timezaps[0] = malloc((arr_size)* sizeof(int));
	timezaps[1] = malloc((arr_size)* sizeof(int));



	if (killswitch){
	  printf("Using file %s as sample killer...\n",killfile);
	  FILE *killer;
	  if ((killer = fopen(killfile, "r")) == NULL){
	    printf("Failed to open the killfile: %s\n",killfile);
	    exit(-1);
	  }
	  int a, j=0;
	  //	  while (!feof(killer)){
	  //	    fscanf(killer,"%i",&a);
	  //	    j++;
	  //	  }
	  //	  fclose(killer);

	  
	  if ((killer = fopen(killfile, "r")) == NULL){
	    printf("Failed to open the killfile: %s\n",killfile);
	    exit(-1);
	  }
	  j=0;
	  char byte;
	  /* Skip the first line of the file */
	  while (!feof(killer)){
	    fread(&byte,1,1,killer);
	    
	    if( byte == '\n' )break;
	  }
	  int good_sample=1;
	  while (!feof(killer)){
	    
	    fscanf(killer,"%i\n",&good_sample);
	    if(!good_sample){
	      if ( ntimezaps > arr_size){
		arr_size*=2;
		timezaps[0]=realloc(timezaps[0],arr_size*sizeof(int));
		timezaps[1]=realloc(timezaps[1],arr_size*sizeof(int));
	      }
	      timezaps[0][ntimezaps]=j;
	      timezaps[1][ntimezaps]=j+1;
	      //	      printf("%d\t%d\t%d\t%d\n",ntimezaps,j,timezaps[0][0],timezaps[1][0]);
	      ntimezaps++;
	    }
	    j++;
	  }
	  fclose(killer);
	  printf("Zapping %d/%d time samples\n",ntimezaps,j);
	}
	
	
	if (optind >= argc) {
		printf("Error: Need a fil file to work on... (use --help for options)\n");
		exit(1);
	}

	//printf("%s\n",argv[optind]);
	if ((file = fopen(argv[optind],"r+")) ==NULL){
	  printf("Failed to open file passed as CLA.\n");
	  exit(-5);
	}

	if ( newname[0]!='\0' || newra < 900000000 || newdec < 900000000)fix_header(file,newname,newra,newdec);

	if(ntimezaps > 0 || nfreqzaps > 0)zap_em(file,timezaps,ntimezaps,fzaps,nfreqzaps,mean,sigma);

	return 0;
}


void fix_header(FILE* file, char* newname, double newra, double newdec){
	int newlen;
	int hdr_len;
	char* hdr_arr;
	char* ptr;
	char buf[1024];
	int an_int;
	double a_double;
	float a_float;
	int i;

	printf("Fixing header\n");
	newlen = strlen(newname);

	fseek(file,0,SEEK_SET);

	hdr_len=read_header(file);
	fseek(file,0,SEEK_SET);

	hdr_arr = (char*)malloc(hdr_len);

	fread(hdr_arr,1,hdr_len,file);

	ptr = hdr_arr;

	while((ptr-hdr_arr) < hdr_len){
		memcpy(buf,ptr,11);
		buf[11]='\0';
		if(newname[0]!='\0' && strcmp(buf,"source_name")==0){
			ptr+=11;
			an_int = *((int*)(ptr));
			ptr+=sizeof(int);
			memcpy(buf,ptr,an_int);
			buf[an_int]='\0';
			printf("old src name = '%s'\n",buf);
			if(an_int > newlen){
				// the old name is longer than the new, pad
				memcpy(buf,newname,newlen);
				for (i=newlen; i < an_int; i++){
					buf[i]=' ';
				}
			} else {
				memcpy(buf,newname,an_int);
			}
			buf[an_int] = '\0';
			printf("new src name = '%s'\n",buf);
			memcpy(ptr,buf,an_int);


		}

		memcpy(buf,ptr,7);
		buf[7]='\0';
		if(newra < 900000000 &&  strcmp(buf,"src_raj")==0){
			ptr+=7;
			a_double = *((double*)(ptr));
			printf("old ra = '%lf'\n",a_double);
			printf("new ra = '%lf'\n",newra);
			*((double*)(ptr)) = newra;
		}

		if(newdec < 900000000 &&  strcmp(buf,"src_dej")==0){
			ptr+=7;
			a_double = *((double*)(ptr));
			printf("old dec = '%lf'\n",a_double);
			printf("new dec = '%lf'\n",newdec);
			*((double*)(ptr)) = newdec;
		}
		ptr++;
	}

	// now re-write the header
	

	fseek(file,0,SEEK_SET);
	fwrite(hdr_arr,1,hdr_len,file);

	free(hdr_arr);
}

#define ARRL 100000
void zap_em(FILE* file, int* tzaps[2], int ntzaps, int fzaps[1024], int nfzaps,float mean, float sigma){
	long long unsigned cur_sample;
	int sample_nbytes;
	int cur_tzap;
	int tz_start;
	int tz_end;
	int c,i;
	int rval;
	char byte;
	int max, min;
	char rndarr[ARRL];
	int mask;
	int rem,stor;

	printf("Zapping data\n");

	byte=0;
	// rewind the file
	fseek(file,0,SEEK_SET);

	// read the file header
	read_header(file);

	if ((nchans*nbits)%8){
		fprintf(stderr,"ERROR: bytes per sample is not an integer\n");
		exit(1);
	}
	sample_nbytes = (nchans*nbits)/8;
	max=(int)(pow(2,nbits))-1;
	min=0;

	cur_sample=0;
	cur_tzap=0;

	if (cur_tzap < ntzaps){
		tz_start = tzaps[0][0];
		tz_end = tzaps[1][0];
	}

	srand ( time(NULL) );
	rem=0;
	for(c=0;c<ARRL;c++){
		byte=0;
		i=0;
		while(1){
			if(i==8){
				rndarr[c] = byte;
				break;
			}

			rval=(int)rint(get_random_value(mean,sigma));
			stor=rval;
			rval+=rem;
			if(rval > max){
				rem = rval-max;
				rval=max;
			} else if(rval < min){
				rem = rval-min;
				rval=min;
			} else {
				rem=0;
			}
			//printf("%d %d %d\n",rval,stor,rem);
			byte |= (rval << i);
			i+=nbits;
		}
	}

	while (cur_tzap < ntzaps){
		while (cur_sample < tz_start){
			// need to fzap, but for now we seek over...
			fseek(file,sample_nbytes,SEEK_CUR);
			cur_sample++;
		}
//		printf("Zapping: %d -> %d\n",tz_start,tz_end);
		// start zapping
		while (cur_sample < tz_end){
		  //		  printf("cur = %d\t%d\t%d\n",cur_sample,tz_start, tz_end);
			i=0;
			for(c=0;c<sample_nbytes;c++){
				i=(int)((rand()*(ARRL-1.0))/(float)RAND_MAX);
				byte=rndarr[i];
				fwrite(&byte,1,1,file);
			}
			cur_sample++;
		}
		cur_tzap++;
		if (cur_tzap < ntzaps){
			tz_start = tzaps[0][cur_tzap];
			tz_end = tzaps[1][cur_tzap];
		}
	}
}

float get_random_value(float mu, float sigma){
	float rnd1,rnd2,ret;
	rnd1 = ((float)rand())/(float)RAND_MAX;
	rnd2 = ((float)rand())/(float)RAND_MAX;
	if(rnd1==0.0) rnd1 = ((float)rand())/(float)RAND_MAX;
	ret= ((sigma * sqrt(- (2.0 * log(rnd1)))) * cos((2.0 * 3.14159265) * rnd2)) + mu;
	return ret;
}
