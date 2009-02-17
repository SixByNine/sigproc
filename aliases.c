#include <string.h>

/* returns the correct TEMPO site code given SIGPROC's telescope_id */
char tempo_site(int telescope_id) /*includefile*/
{
  switch (telescope_id) {
  case 0:
    return('3');
    break;
  case 1:
    return('3'); /*AO*/
    break;
  case 3:
    return('f'); /*Nancay*/
    break;
  case 4:
    return('7'); /*Parkes*/
    break;
  case 5:
    return('8'); /*Jodrell*/
    break;
  case 6:
    return('1'); /*GBT*/
    break;
  case 8:
    return('g'); /*Effelsberg*/
    break;
  default:
    return('?'); /*unknown*/
    /*error_message("tempo_site: unknown telescope!");*/
    break;
  }
}

char *telescope_name (int telescope_id) /*includefile*/
{
  char *telescope,string[80];
  switch (telescope_id) {
  case 0: 
    strcpy(string,"Fake");
    break;
  case 1: 
    strcpy(string,"Arecibo");
    break;
  case 2: 
    strcpy(string,"Ooty");
    break;
  case 3: 
    strcpy(string,"Nancay");
    break;
  case 4: 
    strcpy(string,"Parkes");
    break;
  case 5: 
    strcpy(string,"Jodrell");
    break;
  case 6: 
    strcpy(string,"GBT");
    break;
  case 7: 
    strcpy(string,"GMRT");
    break;
  case 8: 
    strcpy(string,"Effelsberg");
    break;
  default: 
    strcpy(string,"???????"); 
    break;
  }
  telescope=(char *) malloc(strlen(string));
  strcpy(telescope,string);
  return(telescope);
}
char *backend_name (int machine_id) /*includefile*/
{
  char *backend, string[80];
  switch (machine_id) {
  case 0:
    strcpy(string,"FAKE");
    break;
  case 1:
    strcpy(string,"PSPM");
    break;
  case 2:
    strcpy(string,"WAPP");
    break;
  case 3:
    strcpy(string,"AOFTM");
    break;
  case 4:
    strcpy(string,"BPP");
    break;
  case 5:
    strcpy(string,"OOTY");
    break;
  case 6:
    strcpy(string,"SCAMP");
    break;
  case 7:
    strcpy(string,"GMRTFB");
    break;
  case 8:
    strcpy(string,"PULSAR2000");
    break;
  default:
    strcpy(string,"?????");
    break;
  }
  backend=(char *) malloc(strlen(string));
  strcpy(backend,string);
  return(backend);
}
char *data_category (int data_type) /*includefile*/
{
  char *datatype, string[80];
  switch (data_type) {
  case 0:
    strcpy(string,"raw data");
    break;
  case 1:
    strcpy(string,"filterbank");
    break;
  case 2:
    strcpy(string,"time series");
    break;
  case 3:
    strcpy(string,"pulse profiles");
    break;
  case 4:
    strcpy(string,"amplitude spectrum");
    break;
  case 5:
    strcpy(string,"complex spectrum");
    break;
  case 6:
    strcpy(string,"dedispersed subbands");
    break;
  default:
    strcpy(string,"unknown!");
    break;
  }
  datatype=(char *) malloc(strlen(string));
  strcpy(datatype,string);
  return(datatype);
}
