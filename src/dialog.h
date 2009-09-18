
//
// dialog.h (quick replacement for dialog)
// 

#define NELEMENTS 100

class button {

 public:
  float x;
  float y;
  float xmin, xmax;
  float ymin, ymax;
  char * label;
  int pressed;
  int inside(float x, float y);
  int draw();
  button();
  button(float xx, float yy, char * input);
};

class check {
 public:
  float x;
  float y;
  float xmin,xmax,ymin,ymax;
  char * label;
  int on;
  int inside(float x, float y);
  int draw();
  int toggle();
  check();
  check(float xx, float yy, char * input, int status);
};

class radio {
 public:
  float x;
  float y;
  char * label;
  int on;
  int groupid;
  int draw();
  int inside(float x, float y);
  radio();
  radio(float xx, float yy, char * input, int status, int group_input);
};

class plotregion {
 public:
  float xmin;
  float xmax;
  float ymin;
  float ymax;

  float xminworld;
  float xmaxworld;
  float yminworld;
  float ymaxworld;

  int xlog;
  int ylog;

  //   These are set up by cmbpt
  int nplotted;  // number of points in the plot
  float * xdat;  //  x array of plotted data
  float * ydat;  //  y array of plotted data
  int * index;     // index back to parent 

  int active;
  int inside(float *x, float *y);
  int worldcoords(float * x, float * y);
  int set(float xmin, float xmax, float ymin, float ymax);
  int reset();
  int query(float * xmn, float * xmx, float * ymn, float * ymx);
  int erase();
  int erase(float,float,float,float);
};

class staticText{
 public:
  float x;
  float y;
  float orientation;
  char * label;
  int ci;
  float size;
  int draw();
  staticText();
  staticText(float xx, float yy, char * input, int colour, float ss,
	     float o=0.0);
};

class dialog {

 public:
  char * label;
  int nbutton;
  int nradio;
  int ncheck;
  int nplotregion;
  int nstaticText;
  int addbutton(float x, float y, char * label);
  int addradio(float x, float y, char * label, int c, int s);
  int addcheck(float x, float y, char * label, int c);
  int addplotregion(float xmin, float xmax, float ymin, float ymax);
  int addstaticText(float x, float y, char * l, int ci, float size, float o);
  int groupon(int group);
  button * buttons;
  check * checks;
  radio * radios;
  plotregion * plotregions;
  staticText * staticTexts;
  dialog();
  int draw();
  int update();
  int manage(float * x, float * y, char * ans, int * plotno); // returns button id
};

