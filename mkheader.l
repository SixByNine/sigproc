%{
#include "y.tab.h"
#define SKIP(x) if( cf ) return( x )
extern int linecount;
int cf = 1;
double atof();
%}
%%
struct                    { SKIP(STRUCT);   }
int                       { yylval.tok = INTEGER; SKIP(INTEGER);  }
long                      { yylval.tok = LONG;    SKIP(LONG);     }
long[ \t]+long            { yylval.tok = LONGLONG;  SKIP(LONGLONG);     }
double                    { yylval.tok = DOUBLE;  SKIP(DOUBLE);   }
float                     { yylval.tok = FLOAT;   SKIP(FLOAT);    }
char                      { yylval.tok = CHARSTAR; SKIP(CHARSTAR); }
short                     { SKIP(SHORT); }
unsigned                  { SKIP(UNSIGNED); }
\;                        { SKIP(SEMI); }
\{                        { SKIP(LBRACE); }
\}                        { SKIP(RBRACE); }
\[                        { SKIP(LB); }
\]                        { SKIP(RB); }
\/\*[^\n]+\*\/                { yylval.var = (char *)malloc( strlen(yytext)+1);
                            strcpy( yylval.var, yytext);cf=1;return( COMMENT );}
\/\*[^*]*\n               { cf = 0; linecount++; }
\*\/                      { cf = 1; }
extern.*\n                { linecount++; /* ignore */ }
[ \t]                     { ; /* ignore */ }
[A-Za-z0-9_*]+ { yylval.var = (char *)malloc( strlen(yytext)+1);
                       strcpy( yylval.var, yytext );SKIP( VAR );}
\#.*\n              { linecount++; /* ignore */ }
[\n]               { linecount++; ignore(); /* ignore */ }
[^ ]              { SKIP( yytext[0] ); }
%%
