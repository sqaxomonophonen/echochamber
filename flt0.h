#ifndef FLT0_H

#define FLT0_BZ0 (100)
#define FLT0_BEGIN (101)
#define FLT0_END (102)
#define FLT0_SUBSEP (103)
#define FLT0_SEP (104)
#define FLT0_FLOAT (105)

struct flt0_bz0_point {
	float frq,att,pha;
};

static inline const char* flt0_token_cstr(int tok)
{
	switch (tok) {
		case FLT0_BZ0: return "FLT0_BZ0";
		case FLT0_BEGIN: return "FLT0_BEGIN";
		case FLT0_END: return "FLT0_END";
		case FLT0_SUBSEP: return "FLT0_SUBSEP";
		case FLT0_SEP: return "FLT0_SEP";
		case FLT0_FLOAT: return "FLT0_FLOAT";
	}
	return "?";
}

static inline void flt0_unexpected_token(int tok)
{
	fprintf(stderr, "flt0 parse error; unexpected token %s\n", flt0_token_cstr(tok));
	exit(EXIT_FAILURE);
}

float flt0_get_floatval();

#define FLT0_H
#endif
