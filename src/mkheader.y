%{

/*
 * Cactus File @(#)mkgs.y	1.2
 *         SID 1.2
 *        Date 8/6/96
 */

#include <stdio.h>
#include <ctype.h>

%}

%union {
         char *var;
         int tok;
       }

%start list
%token STRUCT LBRACE RBRACE SEMI LB RB 
%token <tok> DOUBLE INTEGER LONG LONGLONG FLOAT SHORT UNSIGNED CHARSTAR BYTE
%token <var> VAR COMMENT

%type <var> construct
%type <tok> notchar
%type <tok> typepp

%%

list :
     | list construct
       { got_construct($2); }
     | list error 
        { yyerrok; yyclearin; YYACCEPT; }
     | list COMMENT 
     ;

construct : STRUCT VAR LBRACE vars RBRACE SEMI
       { $$ = $2; }
      ;

vars : vars line 
     | line
     ;

line : COMMENT
       { comment($1); }
     | CHARSTAR VAR LB VAR RB SEMI
       { add_char($2, $4 ); }
     | notchar VAR SEMI
       { add_notchar( $1, $2); }
     | notchar VAR LB VAR RB SEMI
       { add_array($1, $2, $4); }
     | notchar VAR LB VAR RB LB VAR RB SEMI
       { printf( "/* skipping %s */\n", $2 );}
     | STRUCT VAR VAR LB VAR RB SEMI
       { printf( "/* skipping struct array %s %s */\n", $2, $3 );}
     | STRUCT VAR VAR SEMI
     ;

notchar: typepp
         { $$ = $1; };
       | UNSIGNED typepp
         { $$ = $2; };
       | UNSIGNED CHARSTAR
         { $$ = $2; };
       ;

typepp: DOUBLE
         { $$ = $1; };
       | FLOAT
         { $$ = $1; };
       | LONG
         { $$ = $1; };
       | LONGLONG
         { $$ = $1; };
       | INTEGER
         { $$ = $1; };
       | SHORT
         { $$ = $1; };
       ;


%%

yyerror(s)
char *s;
{
  extern int linecount;

  fprintf( stderr, "%s at line %d\n", s, linecount );
}

yywrap() { return(1); }
