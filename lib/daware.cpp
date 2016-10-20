#include "divsufsort_private.h"
#include "sort/daware.h"

// Call the library
void daware(saidx_t* SAf, saidx_t* SAl, saidx_t* ISAf) {
  sort::suffix::daware(SAf, SAl, ISAf);
}
