/*
**	eraseDP()
**
**	28/08/2002	bklein@MPIfR-Bonn.MPG.de
*/

#include <ctype.h>

void eraseDP( char *cbuf )	/* includefile */
{
  unsigned int	i = 0;

  while( cbuf[i] != '\0' )  {
    if (cbuf[i] == ':')	   cbuf[i] = ' ';
    if (iscntrl(cbuf[i]))  cbuf[i] = '\0';
    i++;
  }
}


