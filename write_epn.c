#include <stdio.h>
#include "epn.h"
/* write an epn file */
void write_epn(FILE *fptr, struct EPN epn) /*includefile*/
{
  int counter,nlines,i;

  nlines=epn.nbins/20;
  if ((epn.nbins-epn.nbins/20)>0) nlines++;
  counter=6+epn.npol*(nlines+2);
  fprintf(fptr,"EPN 6.00%4d%68s",counter,epn.history);
  fprintf(fptr,"%12s%12s%16.12f%8.3f%10.3f%6s%8s%8s",
  epn.jname,epn.cname,epn.pbar,epn.dm,epn.rm,epn.catref,epn.bibref," ");

  fprintf(fptr,"%10.3f%11.3f%-8s%10.3f%8.3f%1c%1c%31s",
  epn.raj,epn.dec,epn.telname,epn.epoch,epn.opos,epn.paflag,epn.timflag," "); 

  fprintf(fptr,"%17.5f%17.5f%17.5f%29s",epn.xtel,epn.ytel,epn.ztel," ");

  fprintf(fptr,"%02d%02d%4d%04d%04d%02d%04d%04d%12.6f%12.6f%06d%04d%4d%1c%15s",
  epn.day,epn.month,epn.year,epn.scanno,epn.subscan,epn.npol,epn.nfreq,
  epn.nbins,epn.tbin,epn.tres,epn.nint,epn.ncal,epn.lcal,epn.fluxflag, " ");
  
  fprintf(fptr,"--------------------------------------------------------------------------------");

  fprintf(fptr,"%8s%4d%4d%12.6f%8s%12.6f%8s%17.5f%7s",
  epn.idfield,epn.nband,epn.navg,epn.f0,epn.uf,epn.df,epn.ud,epn.tstart," ");

  fprintf(fptr,"%#12G%#12G%#12G%16.12f%28s",
  epn.scale,epn.offset,epn.rms,epn.papp," ");

  for (i=0;i<epn.nbins;i++) fprintf(fptr,"%04X",epn.iprofile[i]);
  for (i=1;i<=20-(epn.nbins-(epn.nbins/20)*20);i++) fprintf(fptr,"0000");
}
/* write an epn header */
void write_epn_header(FILE *fptr, struct EPN epn) /*includefile*/
{
  int counter,nlines,i;

  nlines=epn.nbins/20;
  if ((epn.nbins-epn.nbins/20)>0) nlines++;
  counter=6+epn.npol*epn.nfreq*(nlines+2);
  fprintf(fptr,"EPN 6.00%4d%68s",counter,epn.history);
  fprintf(fptr,"%12s%12s%16.12f%8.3f%10.3f%6s%8s%8s",
  epn.jname,epn.cname,epn.pbar,epn.dm,epn.rm,epn.catref,epn.bibref," ");

  fprintf(fptr,"%10.3f%11.3f%-8s%10.3f%8.3f%1c%1c%31s",
  epn.raj,epn.dec,epn.telname,epn.epoch,epn.opos,epn.paflag,epn.timflag," "); 

  fprintf(fptr,"%17.5f%17.5f%17.5f%29s",epn.xtel,epn.ytel,epn.ztel," ");

  fprintf(fptr,"%02d%02d%4d%04d%04d%02d%04d%04d%12.6f%12.6f%06d%04d%4d%1c%15s",
  epn.day,epn.month,epn.year,epn.scanno,epn.subscan,epn.npol,epn.nfreq,
  epn.nbins,epn.tbin,epn.tres,epn.nint,epn.ncal,epn.lcal,epn.fluxflag, " ");
  fprintf(fptr,"--------------------------------------------------------------------------------");
}
/* write an epn file */
void write_epn_subheader(FILE *fptr, struct EPN epn) /*includefile*/
{
  int counter,nlines,i;

  nlines=epn.nbins/20;
  if ((epn.nbins-epn.nbins/20)>0) nlines++;
  counter=6+epn.npol*(nlines+2);

  fprintf(fptr,"%8s%4d%4d%12.6f%8s%12.6f%8s%17.5f%7s",
  epn.idfield,epn.nband,epn.navg,epn.f0,epn.uf,epn.df,epn.ud,epn.tstart," ");

  fprintf(fptr,"%#12G%#12G%#12G%16.12f%28s",
  epn.scale,epn.offset,epn.rms,epn.papp," ");

  for (i=0;i<epn.nbins;i++) fprintf(fptr,"%04X",epn.iprofile[i]);
  for (i=1;i<=20-(epn.nbins-(epn.nbins/20)*20);i++) fprintf(fptr,"0000");
}
