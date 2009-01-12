#include <stdio.h>
int file_exists(char *filename) /*includefile*/
{
  if ((fopen(filename,"rb"))==NULL) {
	return(0);
  } else {
	return(1);
  }
}
