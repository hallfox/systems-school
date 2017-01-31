#include "abstract_matrix.h"
#include "dense_matrix.h"

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>

//TODO: Add types, data and functions as required.

/** Return a newly allocated matrix with all entries in consecutive
 *  memory locations (row-major layout).  All entries in the newly
 *  created matrix are initialized to 0.  Set *err to EINVAL if nRows
 *  or nCols <= 0, to ENOMEM if not enough memory.
 */
DenseMatrix *
newDenseMatrix(int nRows, int nCols, int *err)
{
  return NULL;
}

/** Return implementation of functions for a dense matrix; these functions
 *  can be used by sub-classes to inherit behavior from this class.
 */
const DenseMatrixFns *
getDenseMatrixFns(void)
{
  return NULL;
}
