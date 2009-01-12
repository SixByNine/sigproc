/* ooty header file information */

#define BUFSZ		(256*1024)
#define TAPE_BLOCKS	8
#define	TAPERECSZ	(BUFSZ/TAPE_BLOCKS)

#define	HDRSZ		2048
#define DATASZ		(BUFSZ - HDRSZ)
#define NOBUFS		11

struct dBuf {
  int	seq, owner, id, t_first,t_last ;
  char    project[32];
  char    spare[12];
  char    hdr[HDRSZ-64];
  char    data[DATASZ];
} ;
        

