#ifndef lint
static char const 
yyrcsid[] = "$FreeBSD: src/usr.bin/yacc/skeleton.c,v 1.28 2000/01/17 02:04:06 bde Exp $";
#endif
#include <stdlib.h>
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYLEX yylex()
#define YYEMPTY -1
#define yyclearin (yychar=(YYEMPTY))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING() (yyerrflag!=0)
static int yygrowstack();
#define YYPREFIX "yy"
#line 2 "mkheader.y"

/*
 * Cactus File @(#)mkgs.y	1.2
 *         SID 1.2
 *        Date 8/6/96
 */

#include <stdio.h>
#include <ctype.h>

#line 14 "mkheader.y"
typedef union {
         char *var;
         int tok;
       } YYSTYPE;
#line 33 "y.tab.c"
#define YYERRCODE 256
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
const short yylhs[] = {                                        -1,
    0,    0,    0,    0,    1,    4,    4,    5,    5,    5,
    5,    5,    5,    5,    2,    2,    2,    3,    3,    3,
    3,    3,    3,
};
const short yylen[] = {                                         2,
    0,    2,    2,    2,    6,    2,    1,    1,    6,    3,
    6,    9,    7,    4,    1,    2,    2,    1,    1,    1,
    1,    1,    1,
};
const short yydefred[] = {                                      1,
    0,    3,    0,    4,    2,    0,    0,    0,   18,   22,
   20,   21,   19,   23,    0,    0,    8,    0,   15,    0,
    7,    0,   17,   16,    0,    0,    0,    6,    0,    0,
   10,    0,    5,   14,    0,    0,    0,    0,    0,    0,
    0,    9,   11,    0,   13,    0,    0,   12,
};
const short yydgoto[] = {                                       1,
    5,   18,   19,   20,   21,
};
const short yysindex[] = {                                      0,
 -256,    0, -269,    0,    0, -253, -241, -242,    0,    0,
    0,    0,    0,    0, -230, -233,    0, -229,    0, -255,
    0, -228,    0,    0, -216, -254, -214,    0, -240, -225,
    0, -224,    0,    0, -223, -243, -231, -212, -209, -219,
 -208,    0,    0, -218,    0, -207, -204,    0,
};
const short yyrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,
};
const short yygindex[] = {                                      0,
    0,    0,   38,    0,   37,
};
#define YYTABLESIZE 57
const short yytable[] = {                                       2,
    3,    8,    6,   27,    7,   31,   32,    9,   10,   11,
   12,   13,   14,   15,   16,    8,    4,   17,   39,   34,
   35,    9,   10,   11,   12,   13,   14,   15,   16,   22,
   40,   17,    9,   10,   11,   12,   13,   14,   25,   23,
   43,   44,   26,   29,   30,   33,   36,   37,   38,   41,
   42,   45,   24,   46,   47,   48,   28,
};
const short yycheck[] = {                                     256,
  257,  257,  272,  259,  258,  260,  261,  263,  264,  265,
  266,  267,  268,  269,  270,  257,  273,  273,  262,  260,
  261,  263,  264,  265,  266,  267,  268,  269,  270,  272,
  262,  273,  263,  264,  265,  266,  267,  268,  272,  270,
  260,  261,  272,  272,  261,  260,  272,  272,  272,  262,
  260,  260,   15,  272,  262,  260,   20,
};
#define YYFINAL 1
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 273
#if YYDEBUG
const char * const yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"STRUCT","LBRACE","RBRACE","SEMI",
"LB","RB","DOUBLE","INTEGER","LONG","LONGLONG","FLOAT","SHORT","UNSIGNED",
"CHARSTAR","BYTE","VAR","COMMENT",
};
const char * const yyrule[] = {
"$accept : list",
"list :",
"list : list construct",
"list : list error",
"list : list COMMENT",
"construct : STRUCT VAR LBRACE vars RBRACE SEMI",
"vars : vars line",
"vars : line",
"line : COMMENT",
"line : CHARSTAR VAR LB VAR RB SEMI",
"line : notchar VAR SEMI",
"line : notchar VAR LB VAR RB SEMI",
"line : notchar VAR LB VAR RB LB VAR RB SEMI",
"line : STRUCT VAR VAR LB VAR RB SEMI",
"line : STRUCT VAR VAR SEMI",
"notchar : typepp",
"notchar : UNSIGNED typepp",
"notchar : UNSIGNED CHARSTAR",
"typepp : DOUBLE",
"typepp : FLOAT",
"typepp : LONG",
"typepp : LONGLONG",
"typepp : INTEGER",
"typepp : SHORT",
};
#endif
#if YYDEBUG
#include <stdio.h>
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH 10000
#endif
#endif
#define YYINITSTACKSIZE 200
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short *yyss;
short *yysslim;
YYSTYPE *yyvs;
int yystacksize;
#line 85 "mkheader.y"

yyerror(s)
char *s;
{
  extern int linecount;

  fprintf( stderr, "%s at line %d\n", s, linecount );
}

yywrap() { return(1); }
#line 188 "y.tab.c"
/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack()
{
    int newsize, i;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = yystacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;
    i = yyssp - yyss;
    newss = yyss ? (short *)realloc(yyss, newsize * sizeof *newss) :
      (short *)malloc(newsize * sizeof *newss);
    if (newss == NULL)
        return -1;
    yyss = newss;
    yyssp = newss + i;
    newvs = yyvs ? (YYSTYPE *)realloc(yyvs, newsize * sizeof *newvs) :
      (YYSTYPE *)malloc(newsize * sizeof *newvs);
    if (newvs == NULL)
        return -1;
    yyvs = newvs;
    yyvsp = newvs + i;
    yystacksize = newsize;
    yysslim = yyss + newsize - 1;
    return 0;
}

#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab

#ifndef YYPARSE_PARAM
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG void
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif	/* ANSI-C/C++ */
#else	/* YYPARSE_PARAM */
#ifndef YYPARSE_PARAM_TYPE
#define YYPARSE_PARAM_TYPE void *
#endif
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG YYPARSE_PARAM_TYPE YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL YYPARSE_PARAM_TYPE YYPARSE_PARAM;
#endif	/* ANSI-C/C++ */
#endif	/* ! YYPARSE_PARAM */

int
yyparse (YYPARSE_PARAM_ARG)
    YYPARSE_PARAM_DECL
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register const char *yys;

    if ((yys = getenv("YYDEBUG")))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    if (yyss == NULL && yygrowstack()) goto yyoverflow;
    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if ((yyn = yydefred[yystate])) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yysslim && yygrowstack())
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#if defined(lint) || defined(__GNUC__)
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#if defined(lint) || defined(__GNUC__)
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yysslim && yygrowstack())
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 2:
#line 32 "mkheader.y"
{ got_construct(yyvsp[0].var); }
break;
case 3:
#line 34 "mkheader.y"
{ yyerrok; yyclearin; YYACCEPT; }
break;
case 5:
#line 39 "mkheader.y"
{ yyval.var = yyvsp[-4].var; }
break;
case 8:
#line 47 "mkheader.y"
{ comment(yyvsp[0].var); }
break;
case 9:
#line 49 "mkheader.y"
{ add_char(yyvsp[-4].var, yyvsp[-2].var ); }
break;
case 10:
#line 51 "mkheader.y"
{ add_notchar( yyvsp[-2].tok, yyvsp[-1].var); }
break;
case 11:
#line 53 "mkheader.y"
{ add_array(yyvsp[-5].tok, yyvsp[-4].var, yyvsp[-2].var); }
break;
case 12:
#line 55 "mkheader.y"
{ printf( "/* skipping %s */\n", yyvsp[-7].var );}
break;
case 13:
#line 57 "mkheader.y"
{ printf( "/* skipping struct array %s %s */\n", yyvsp[-5].var, yyvsp[-4].var );}
break;
case 15:
#line 62 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 16:
#line 64 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 17:
#line 66 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 18:
#line 70 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 19:
#line 72 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 20:
#line 74 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 21:
#line 76 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 22:
#line 78 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
case 23:
#line 80 "mkheader.y"
{ yyval.tok = yyvsp[0].tok; }
break;
#line 455 "y.tab.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yysslim && yygrowstack())
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
