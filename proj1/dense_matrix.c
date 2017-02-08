#include "abstract_matrix.h"
#include "dense_matrix.h"

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>

typedef struct {
  DenseMatrix;
  int nRows;
  int nCols;
  MatrixBaseType mat[];
} DenseMatrixImpl;

/** Examines the matrix as a DenseMatrix, and verifies that it is
    in a valid state, otherwise, set *err to EINVAL. */
static void verifyDenseMatrix(const Matrix *this, int *err)
{
  const DenseMatrixImpl *matrix = (const DenseMatrixImpl *)this;
  if (matrix->nRows <= 0 || matrix->nCols <= 0) {
    *err = EINVAL;
  }
}

static const char *getKlass(const Matrix *this, int *err)
{
  verifyDenseMatrix(this, err);
  return "denseMatrix";
}

static void freeDenseMatrix(Matrix *this, int *err)
{
  verifyDenseMatrix(this, err);
  DenseMatrixImpl *matrix = (DenseMatrixImpl *)this;
  free(matrix);
}

static int getNRows(const Matrix *this, int *err)
{
  verifyDenseMatrix(this, err);
  const DenseMatrixImpl *matrix = (const DenseMatrixImpl *)this;
  return matrix->nRows;
}

static int getNCols(const Matrix *this, int *err)
{
  verifyDenseMatrix(this, err);
  const DenseMatrixImpl *matrix = (const DenseMatrixImpl *)this;
  return matrix->nCols;
}

static MatrixBaseType getElement(const Matrix *this,
				 int rowIndex, int colIndex, int *err)
{
  const DenseMatrixImpl *matrix = (const DenseMatrixImpl *)this;
  int nRows = this->fns->getNRows(this, err);
  if (*err == EINVAL) return 0;
  int nCols = this->fns->getNCols(this, err);
  if (*err == EINVAL) return 0;
  // Range check
  if (rowIndex >= nRows || colIndex >= nCols) {
    *err = EDOM;
    return 0;
  }

  return matrix->mat[rowIndex*nCols+colIndex];
}

static void setElement(Matrix *this, int rowIndex, int colIndex,
		       MatrixBaseType element, int *err)
{
  DenseMatrixImpl *matrix = (DenseMatrixImpl *)this;
  int nRows = this->fns->getNRows(this, err);
  if (*err == EINVAL) return;
  int nCols = this->fns->getNCols(this, err);
  if (*err == EINVAL) return;
  // Range check
  if (rowIndex >= nRows || colIndex >= nCols) {
    *err = EDOM;
    return;
  }

  matrix->mat[rowIndex*nCols+colIndex] = element;
}

static _Bool isInit = false;
static DenseMatrixFns denseMatrixFns = {
  .getKlass = getKlass,
  .free = freeDenseMatrix,
  .getNRows = getNRows,
  .getNCols = getNCols,
  .getElement = getElement,
  .setElement = setElement,
};

static void patchDenseMatrixFns(void)
{
  // Set up the methods if not already done
  if (!isInit) {
    const MatrixFns *fns = getAbstractMatrixFns();
    denseMatrixFns.transpose = fns->transpose;
    denseMatrixFns.mul = fns->mul;
    isInit = true;
  }
}

/** Return a newly allocated matrix with all entries in consecutive
 *  memory locations (row-major layout).  All entries in the newly
 *  created matrix are initialized to 0.  Set *err to EINVAL if nRows
 *  or nCols <= 0, to ENOMEM if not enough memory.
 */
DenseMatrix *
newDenseMatrix(int nRows, int nCols, int *err)
{
  // Check if dimensions make sense
  if (nRows <=0 || nCols <= 0) {
    *err = EINVAL;
    return NULL;
  }

  // Allocate the matrix with rows, check
  DenseMatrixImpl *matrix = malloc(sizeof(DenseMatrixImpl) +
				   nRows*nCols*sizeof(MatrixBaseType *));
  if (!matrix) {
    *err = ENOMEM;
    return NULL;
  }

  matrix->nRows = nRows;
  matrix->nCols = nCols;
  matrix->fns = (MatrixFns *)getDenseMatrixFns();
  
  return (DenseMatrix *)matrix;
}

/** Return implementation of functions for a dense matrix; these functions
 *  can be used by sub-classes to inherit behavior from this class.
 */
const DenseMatrixFns *
getDenseMatrixFns(void)
{
  patchDenseMatrixFns();
  return &denseMatrixFns;
}
