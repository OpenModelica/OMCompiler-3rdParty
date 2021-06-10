#ifndef DDASKR_TYPES
#define DDASKR_TYPES

/* Table of types */
typedef enum { _FALSE_, _TRUE_ } _boolean;

typedef double real_number;

typedef int integer;
typedef unsigned int uinteger;
#endif

#ifdef __cplusplus
typedef int (*Unknown_fp)(...);
#else
typedef int (*Unknown_fp)();
#endif

/* Table of defines */
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

