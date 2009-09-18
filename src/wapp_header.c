
#include <stdio.h>
#include <fcntl.h>

void head_input( char *, int *, int );
int linecount = 1;

#define YY_INPUT(buf,result,max) head_input(buf,&result,max)
#include "mkheaderlex.c"

#include "key.h"


struct HEADERP *yacc_input = NULL;

/* this routine is not re-entrant, (because of yyparse and yacc_input) */

struct HEADERP *head_parse( FILE *fptr )
{
  int fd, count, ret;
  unsigned char byte;
  struct HEADERKEY *key;
  struct HEADERP *h;

  count = 0;
  while( (ret = fread( &byte, 1, 1, fptr)) == 1 ) {
    if( byte == 0 )
      break;
    count++;
  }

  if( ret < 1 )
    return(NULL);

  h = ( struct HEADERP *)malloc( sizeof(struct HEADERP));
  bzero( h, sizeof(struct HEADERP));

  h->offset = count+1;
  h->fd = fd;
  h->buf = (char *)malloc( h->offset);
  lseek(fd, 0, SEEK_SET);
  read(h->buf, sizeof(h->buf), h->offset, fptr );

  yacc_input = h;

  while(yyparse()); /* use yacc to parse the header */
  yyparse();

  h->headlen = count_size(h);
  h->header = (void *)malloc( h->headlen );
  if( fread( h->header, sizeof(h->header), h->headlen, fptr ) != h->headlen ) {
    perror("read header");
    return(NULL);
  }

  key = h->head;
  count = 0;
  while( key ) {
    key->len = key_sizes(key->type);
    key->offset = count;
    count  += key->len* key->alen;
    key = key->next;
  }

  return(h);
}

count_size(h)
struct HEADERP *h;
{
  struct HEADERKEY *k;
  int count;

  k = h->head;
  count = 0;
  while(k) {
    count += key_sizes(k->type)*k->alen;
    k = k->next;
  }
  return(count);
}

int key_sizes(type)
int type;
{
  switch( type ) {
    case INTEGER:
      return sizeof(int);
    case LONG:
      return sizeof(long);
    case LONGLONG:
      return sizeof(long long);
    case DOUBLE:
      return sizeof(double);
    case FLOAT:
      return sizeof(float);
    case CHARSTAR:
      return sizeof(char);
    default:
      printf("key_sizes error, bad type %d\n", type );
      return 0;
  }
}

int find_hdrval( h, name, hdrval )
struct HEADERP *h;
char *name;
struct HEADERVAL *hdrval;
{
  struct HEADERKEY *key;
  puts("???");
  key = h->head;
  while( key ) {
    if( strcmp( name, key->name ) == 0 )
      break;
    key = key->next;
  }

  if( key ) {
    if( hdrval ) {
      hdrval->value = (void *)((unsigned char *)h->header + key->offset);
      hdrval->key = key;
    }
  } else if( hdrval ) {
    bzero( hdrval, sizeof(struct HEADERVAL));
  }

  return(!key);
}

int set_hdrval( h, name, data, ix )
struct HEADERP *h;
char *name;
void *data;
int ix;
{
  struct HEADERVAL val;
  void *value;

  if( find_hdrval( h, name, &val ) )
    return(-1);

  if( val.key->alen == 1 )
    ix = 0;

  if( ix >= val.key->alen )
    return(-1);

  value = (char *)val.value + val.key->len*ix;

  bcopy( data, value, val.key->len );
  return(0);
}

int write_hdr( h )
struct HEADERP *h;
{
  lseek( h->fd, h->offset, SEEK_SET );
  if( write( h->fd, h->header, h->headlen ) != h->headlen ) {
    perror("write_hdr write");
    return(-1);
  }
  return(0);
}

ignore()
{ 
}


struct HEADERKEY *new_header()
{
  struct HEADERKEY *p;
  extern struct HEADERP *yacc_input;

  p = (struct HEADERKEY *)malloc( sizeof(struct HEADERKEY));
  bzero( p, sizeof(struct HEADERKEY));

  if( !yacc_input->head  )
    yacc_input->head = p;

  if( yacc_input->tail )
    yacc_input->tail->next = p;

  yacc_input->tail = p;
  return(p);
}

add_notchar(v, name)
int v;
char *name;
{
  
  struct HEADERKEY *p = new_header();

  p->name = name;
  p->type = v;
  p->alen = 1;
}


add_array(v, name, sz)
int v;
char *name, *sz;
{
  struct HEADERKEY *p = new_header();

  p->name = name;
  p->type = v;
  if( sz )
    p->alen = atoi(sz);
  else
    p->alen = 1;
}


comment(p)
char *p;
{
  free(p);
}

add_char(name, sz)
char *name, *sz;
{
  struct HEADERKEY *p = new_header();

  p->name = name;
  p->type = CHARSTAR;
  if( sz )
    p->alen = atoi(sz);
  else
    p->alen = 1;
}

got_construct( name )
char *name;
{
  /* printf("got construct\n"); */
}


void head_input( buf, result, max )
char *buf;
int *result;
int max;
{
  extern struct HEADERP *yacc_input;  
  int ct;

  ct = yacc_input->offset -  yacc_input->yacc_offset;
  if( ct <=0 ) {
    *buf = 0;
    *result = 0;
    return;
  }
  if( ct < max ) {
    bcopy( &yacc_input->buf[yacc_input->yacc_offset], buf, ct );
    yacc_input->yacc_offset = yacc_input->offset;
    *result = ct;
  } else { 
    bcopy( &yacc_input->buf[yacc_input->yacc_offset], buf, max );
    yacc_input->yacc_offset += max;
    *result = max;
  }
}

close_parse( h )
struct HEADERP *h;
{
  struct HEADERKEY *key, *fr;

  close(h->fd);

  key = h->head;
  while( key ) {
    fr = key;
    key = key->next;
    free(fr);
  }
  
  free(h->buf);
  free(h->header);
  free(h);
}

int jgetc( FILE *a)
{
  static count = 0; 

  if( count < yacc_input->offset )
    return( yacc_input->buf[count++]);
  else
    return(EOF);
}
