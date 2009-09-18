#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "wapp_header.h"
#include "key.h"

struct WAPP_HEADER head;
/*--------------------------------------------------------------*
 * Generic Conversion from HEADEREP to write into 
 *  a destination.
 *--------------------------------------------------------------*/
int fetch_hdrval(struct HEADERP *h,char *name,void *dest,int ndest) 
{
  struct HEADERVAL hdrval;
  int nsrc;
  bzero(dest,ndest);
  if( find_hdrval(h,name,&hdrval))
  {
    /*fprintf(stderr,"Error Finding header value: %s\n",name);*/
    return(0);
  }
  nsrc = hdrval.key->alen * hdrval.key->len;
  bcopy(hdrval.value,dest,(ndest<nsrc)?ndest:nsrc);
  return(1);
}
