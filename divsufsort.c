#ifdef BUILD_DIVSUFSORT64
#  undef BUILD_DIVSUFSORT64
#endif
#define HAVE_CONFIG_H 1
#include <stdio.h>
#include "divsufsort_private.h"

#include "divsufsort.c.inc.h"
# if (_IS_USED_DAWARE_SORT)
/* daware.cpp */
#   include "daware.cpp.inc.h"
# else
#   include "trsort.c.inc.h"
# endif
#define lg_table sssort_lg_table
#include "sssort.c.inc.h"
#include "utils.c.inc.h"

