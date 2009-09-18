
//
// dialog.C (quick replacement for dialog)
// 

#include "dialog.h"
#include <stdio.h>
#include <stdlib.h>
#include "cpgplot.h"
#include "string.h"
#include "math.h"

int dialog::addbutton(float xx, float yy, char * llabel){

  if (nbutton>NELEMENTS) {
    fprintf(stderr,"Too many buttons in dialog\n");
    exit(-1);
  }
  buttons[nbutton].x = xx;
  buttons[nbutton].y = yy;
  buttons[nbutton].label = new char[strlen(llabel)+1];
  strcpy(buttons[nbutton].label,llabel);
  buttons[nbutton].pressed = 0;
  nbutton++;
  /* return the button number created */
  return(nbutton-1);
}

int dialog::addradio(float xx, float yy, char * llabel, int status, int group){
  
  if (nradio>NELEMENTS) {
    fprintf(stderr,"Too many radios in dialog\n");
    exit(-1);
  }
  radios[nradio].x = xx;
  radios[nradio].y = yy;
  radios[nradio].on = status;
  radios[nradio].groupid = group;
  radios[nradio].label = new char[strlen(llabel)+1];
  strcpy(radios[nradio].label,llabel);
  nradio++;
  return(nradio-1);
}

int dialog::addcheck(float xx, float yy, char * llabel, int status){
  
  if (ncheck>NELEMENTS) {
    fprintf(stderr,"Too many checks in dialog\n");
    exit(-1);
  }
  checks[ncheck].xmin = xx, checks[ncheck].ymin = yy;
  checks[ncheck].xmax = xx+0.02; checks[ncheck].ymax=yy+0.02;
  checks[ncheck].x = xx;
  checks[ncheck].y = yy;
  checks[ncheck].label = new char[strlen(llabel)+1];
  checks[ncheck].on = status;
  strcpy(checks[ncheck].label,llabel);
  ncheck++;
  return(ncheck-1);
}

int dialog::addstaticText(float xx, float yy, char * llabel, int c,
			  float s, float o){
  
  if (nstaticText>NELEMENTS) {
    fprintf(stderr,"Too many staticTexts in dialog\n");
    exit(-1);
  }
  staticTexts[nstaticText].ci=c;
  staticTexts[nstaticText].size=s;
  staticTexts[nstaticText].x = xx;
  staticTexts[nstaticText].y = yy;
  staticTexts[nstaticText].orientation = o;
  staticTexts[nstaticText].label = new char[strlen(llabel)+1];
  strcpy(staticTexts[nstaticText].label,llabel);
  nstaticText++;
  return(0);
}

int dialog::addplotregion(float xmn, float xmx, float ymn, float ymx) {
  
  if (nplotregion>NELEMENTS) {
    fprintf(stderr,"Too many plotregions in dialog\n");
    exit(-1);
  }
  plotregions[nplotregion].xmin = xmn;
  plotregions[nplotregion].xmax = xmx;
  plotregions[nplotregion].ymin = ymn;
  plotregions[nplotregion].ymax = ymx;
  plotregions[nplotregion].xminworld = 0.0;
  plotregions[nplotregion].xmaxworld = 1.0;
  plotregions[nplotregion].yminworld = 0.0;
  plotregions[nplotregion].ymaxworld = 1.0;
  nplotregion++;
  return(0);
}

int dialog::draw(){

  cpgbbuf();
  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0.0,1.0,0.0,1.0);
  for (int i=0;i<nbutton;i++) buttons[i].draw();
  for (int i=0;i<nradio;i++) radios[i].draw();
  for (int i=0;i<ncheck;i++) checks[i].draw();
  for (int i=0;i<nstaticText;i++) staticTexts[i].draw();
  cpgebuf();
  return(0);
}

int dialog::update(){   // Only draws what has changed
  cpgbbuf();
  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0.0,1.0,0.0,1.0);
  for (int i=0;i<nbutton;i++) buttons[i].draw();
  for (int i=0;i<nradio;i++) radios[i].draw();
  for (int i=0;i<ncheck;i++) checks[i].draw();
  for (int i=0;i<nstaticText;i++) staticTexts[i].draw();
  cpgebuf();
  return(0);
}

int dialog::manage(float * x, float * y, char * ans, int * plotno){

  int leave;
  int button = 0;
  *plotno = -1;
  *ans = ' ';
  static float lastx, lasty;

  *x=lastx;
  *y=lasty;

  leave = 0;
  while (!leave){

  // get a keypress
    cpgsvp(0.0,1.0,0.0,1.0);
    cpgswin(0.0,1.0,0.0,1.0);
    cpgcurs(x,y,ans);
    lastx=*x;
    lasty=*y;
    button = -1;
  // First of all, check if it is inside a plot
     for (int i=0;i<nplotregion;i++){
       if (plotregions[i].inside(x,y)){


         button = -1;
         *plotno = i;
	 leave = 1;
        } 
     }
  // Then check if it is inside a button
     for (int i=0;i<nbutton;i++){
       buttons[i].pressed = 0;
       if (buttons[i].inside(*x,*y)){
         button = i;
         buttons[i].pressed = 1;
	 leave = 1;
        } 
     }

  // Then radio controls
     for (int i=0;i<nradio;i++){
       if (radios[i].inside(*x,*y)){
        //  zero all in the group
        for (int j=0;j<nradio;j++) 
          if (radios[j].groupid==radios[i].groupid) {
	    radios[j].on=0;
	    radios[j].draw();
	  }
        // draw this one
	radios[i].on=1;
	  radios[i].draw();
	  leave = 0;
        } 
      }

  // Then check boxes controls
     for (int i=0;i<ncheck;i++){
       if (checks[i].inside(*x,*y)){
         checks[i].toggle();
         checks[i].draw();
	 leave=0;
        } 
     }
  // else stay in loop unless q pressed.
  }
  return(button);
}

int button::inside(float xi, float yi){
if (xi>=xmin && xi<=xmax && yi>=ymin && yi<=ymax) return (1);
     return(0);
}

int button::draw(){
cpgsvp(0.0,1.0,0.0,1.0);
cpgswin(0.0,1.0,0.0,1.0);
cpgsfs(2);

float xl, yl;
cpglen(4,label,&xl,&yl);
xmin = x-2.0*0.005; 
xmax = x + xl + 2.0 * 0.005;
ymin = y-2.0*0.005;
ymax = y+0.015 + 2.0 * 0.005;

cpgsci(1);
cpgrect(x-0.005, x+xl+0.005, y-0.005, y+0.015 + 0.005);
cpgrect(x-2.0*0.005, x+xl+2.0*0.005, y-2.0*0.005, y+0.015 + 2.0*0.005);
cpgtext(x,y,label);
 return(0);
}

button::button(){
  label=NULL;
}

button::button(float xi, float yi, char * l){
x=xi;
y=yi;
label = new char[strlen(l)+1];
strcpy(label,l);
float xl, yl;
cpglen(4,label,&xl,&yl);
xmin = x-2.0*0.005; 
xmax = x + xl + 2.0 * 0.005;
ymin = y-2.0*0.005;
ymax = y+0.015 + 2.0 * 0.005;

}

check::check(){ x=y=0.0;label=NULL; on=0;}

check::check(float xx, float yy, char * l, int status){
  x=xx; y=yy; on=status;
  xmin = xx, ymin = yy;
  xmax = xx+0.02; ymax=yy+0.02;
  label = new char[strlen(l)+1];
  strcpy(label,l);
}

int check::inside(float xx, float yy){
  if (xx>=xmin && xx <=xmax && yy>= ymin && yy<=ymax) return (1);
  return(0);
}

int check::toggle(){
  if (on) {
    on = 0;
    draw();
  }
  else{
    on = 1;
    draw();
  }
  return(0);
}

int check::draw(){

  if (on) {
    cpgsci(2);
    cpgsfs(1);
    cpgrect(xmin,xmax,ymin,ymax);
    cpgsci(1);
    cpgsfs(2);
    cpgrect(xmin,xmax,ymin,ymax);
  }
  else{
    cpgsci(0);
    cpgsfs(1);
    cpgrect(xmin,xmax,ymin,ymax);
    cpgsci(1);
    cpgsfs(2);
    cpgrect(xmin,xmax,ymin,ymax);
  }
  cpgtext(x+0.05,y,label);
  return(0);
}

void fillcircle(float x, float y, int ci){
  float xpts[40], ypts[40];
  for (int i=0;i<40;i++){
    xpts[i]=x+0.01*cos(2.0*M_PI*i/40);
    ypts[i]=y+0.01*sin(2.0*M_PI*i/40);
  }
  cpgsfs(1);
  cpgsci(ci);
  cpgpoly(40,xpts,ypts);
  cpgsfs(2);
  cpgsci(1);
  cpgpoly(40,xpts,ypts);
}

int radio::draw(){  
  if (on) fillcircle(x,y,2); else fillcircle(x,y,0);
  cpgsci(1);
  cpgtext(x+0.02,y-0.005,label);
  return(0);
}

int radio::inside(float xx, float yy){
  if ( (x-xx)*(x-xx) + 
       (y-yy)*(y-yy) <0.01*0.01) {
    return(1);
  }
  return(0);
}

radio::radio(){label=NULL;}
radio::radio(float xx, float yy, char * l, int status, int group_input){
  x = xx; y = yy; label = new char[strlen(l)+1]; strcpy(label,l);
  on=status; groupid = group_input;
}

int plotregion::erase(){
  cpgsfs(1);
  cpgsci(0);
  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0.0,1.0,0.0,1.0);
  cpgrect(xmin-0.065,xmax,ymin-0.1,ymax+0.025);
  cpgsci(1);
  return (0);
}
int plotregion::erase(float dx1, float dx2, float dy1, float dy2){
  cpgsfs(1);
  cpgsci(0);
  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0.0,1.0,0.0,1.0);
  cpgrect(xmin-dx1,xmax+dx2,ymin-dy1,ymax+dy2);
  cpgsci(1);
  return (0);
}

int plotregion::inside(float * x, float * y){
  if (*x>=xmin && *x<=xmax && *y>=ymin && *y<=ymax) {
    *x = (*x-xmin) / ( xmax-xmin) * (xmaxworld-xminworld) + xminworld;
    *y = (*y-ymin) / ( ymax-ymin) * (ymaxworld-yminworld) + yminworld;
    return(1);
  }
  return(0);
}

int plotregion::worldcoords(float * x, float * y){
  *x = (*x-xmin) / ( xmax-xmin) * (xmaxworld-xminworld) + xminworld;
  *y = (*y-ymin) / ( ymax-ymin) * (ymaxworld-yminworld) + yminworld;
  return(0);
}

// Sets the plot coordinates
int plotregion::set(float xmn, float xmx, float ymn, float ymx){
  xminworld = xmn; yminworld = ymn; xmaxworld = xmx; ymaxworld = ymx;
  cpgsvp(xmin,xmax,ymin,ymax);
  cpgswin(xmn,xmx,ymn,ymx);
  return(0);
}

// Sets the plot back to last coordinates
int plotregion::reset(){
  cpgsvp(xmin,xmax,ymin,ymax);
  cpgswin(xminworld,xmaxworld,yminworld,ymaxworld);
  return(0);
}

int plotregion::query(float * xmn, float * xmx, float * ymn, float * ymx){
  *xmn=xminworld;
  *xmx = xmaxworld;
  *ymn = yminworld;
  *ymx = ymaxworld;
  return(0);
}

staticText::staticText(){label=NULL;}
staticText::staticText(float xx, float yy, char * input, int colour, float ss,
		       float o){
  x=xx;y=yy;ci=colour;size=ss; label = new char[strlen(input)+1];
  orientation=o;
  strcpy(label,input);
}

int staticText::draw(){
  cpgsci(ci); cpgsch(size);
  cpgptxt(x,y,orientation,0.0,label);
  cpgsci(1);
  cpgsch(1.0);
  return (0);
}

dialog::dialog(){
  label=NULL;
  nbutton = nradio = ncheck = nstaticText = nplotregion = 0;
  buttons = new button[NELEMENTS];
  radios = new radio[NELEMENTS];
  plotregions = new plotregion[NELEMENTS];
  checks = new check[NELEMENTS];
  staticTexts = new staticText[NELEMENTS];
}

int dialog::groupon(int groupid){
  int elementno;
  elementno=-1;
  for (int i=0;i<nradio;i++) { 
    if (radios[i].groupid==groupid) {
      elementno++;
      if (radios[i].on==1) return(elementno);
    }
  }
  return(-1);
}

