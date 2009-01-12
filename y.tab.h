#ifndef YYERRCODE
#define YYERRCODE 256
#endif

#define STRUCT 257
#define LBRACE 258
#define RBRACE 259
#define SEMI 260
#define LB 261
#define RB 262
#define DOUBLE 263
#define INTEGER 264
#define LONG 265
#define LONGLONG 266
#define FLOAT 267
#define SHORT 268
#define UNSIGNED 269
#define CHARSTAR 270
#define BYTE 271
#define VAR 272
#define COMMENT 273
typedef union {
         char *var;
         int tok;
       } YYSTYPE;
extern YYSTYPE yylval;
