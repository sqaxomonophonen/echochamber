%{
#include "flt0.h"

static double doubleval;
%}

%option nounput
%option noinput

%%

[ \n\t\r] {
}

"bz0" {
	return FLT0_BZ0;
}

"(" {
	return FLT0_BEGIN;
}

")" {
	return FLT0_END;
}

":" {
	return FLT0_SUBSEP;
}

"," {
	return FLT0_SEP;
}

-?([0-9]+|[0-9]*\.[0-9]+) {
	double sign = 1;
	int of = 0;

	if (yytext[0] == '-') {
		sign = -1;
		of++;
	}
	int beg = of;
	while (yytext[of] != 0 && yytext[of] != '.') of++;

	doubleval = 0;

	double m = 1;
	for (int i = of-1; i >= beg; i--) {
		doubleval += (yytext[i] - '0') * m;
		m *= 10;
	}

	if (yytext[of] == '.') {
		of++;
		m = 0.1;
		while (yytext[of] != 0) {
			doubleval += (yytext[of] - '0') * m;
			m *= 0.1;
			of++;
		}
	}

	doubleval *= sign;

	return FLT0_FLOAT;
}

. {
	YY_FATAL_ERROR("flt0 syntax error");
}

%%

int flt0wrap(void)
{
	return 1;
}

float flt0_get_floatval()
{
	return (float) doubleval;
}
