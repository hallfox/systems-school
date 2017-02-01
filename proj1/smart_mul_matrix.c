#include "dense_matrix.h"
#include "smart_mul_matrix.h"

#include <errno.h>
#include <stdbool.h>

typedef struct {
  DenseMatrix;
} SmartMulMatrixImpl;

static const char *getKlass(const Matrix *this, int *err)
{
  return "smartMulMatrix";
}

static void mul(const Matrix *this, const Matrix *multiplier,
		Matrix *product, int *err)
{
  // Check if the dimensions are correct:
  // MxN * NxP = MxP
  const int this_m = this->fns->getNRows(this, err);
  if (*err == EINVAL) return;
  const int this_n = this->fns->getNCols(this, err);
  if (*err == EINVAL) return;
  const int mul_n = multiplier->fns->getNRows(multiplier, err);
  if (*err == EINVAL) return;
  const int mul_p = multiplier->fns->getNCols(multiplier, err);
  if (*err == EINVAL) return;
  const int pr_m = product->fns->getNRows(product, err);
  if (*err == EINVAL) return;
  const int pr_p = product->fns->getNCols(product, err);
  if (*err == EINVAL) return;
  if (!(this_m == pr_m && this_n == mul_n && mul_p == pr_p)) {
    *err = EDOM;
    return;
  }

  // Transpose multiplier, so columns are more likely to end up in the cache
  // NxP -> PxN
  Matrix *tr_multiplier = (Matrix *)newSmartMulMatrix(mul_p, mul_n, err);
  if (*err == EINVAL || *err == ENOMEM) return;
  multiplier->fns->transpose(multiplier, tr_multiplier, err);
  if (*err == EINVAL || *err == EDOM) {
    tr_multiplier->fns->free(tr_multiplier, err);
    return;
  }

  // Do the multiplication
  for (int pr_r = 0; pr_r < pr_m; pr_r++) {
    for (int pr_c = 0; pr_c < pr_p; pr_c++) {
      // Pr[r][c] <- Sum_i This[r][i]*tr_That[c][i]
      MatrixBaseType res = 0;
      for (int i = 0; i < this_n; i++) {
	MatrixBaseType a = this->fns->getElement(this, pr_r, i, err);
	if (*err == EINVAL || *err == EDOM) return;
	MatrixBaseType b = tr_multiplier->fns->getElement(tr_multiplier, pr_c, i, err);
	if (*err == EINVAL || *err == EDOM) return;
	res += a*b;
      }
      product->fns->setElement(product, pr_r, pr_c, res, err);
      if (*err == EINVAL || *err == EDOM) return;
    }
  }

  // Clean up
  tr_multiplier->fns->free(tr_multiplier, err);
}

//TODO: Add types, data and functions as required.
static _Bool isInit = false;
static SmartMulMatrixFns smartMulMatrixFns = {
  .getKlass = getKlass,
  .mul = mul,
};

/** Return a newly allocated matrix with all entries in consecutive
 *  memory locations (row-major layout).  All entries in the newly
 *  created matrix are initialized to 0.  The return'd matrix uses
 *  a smart multiplication algorithm to avoid caching issues;
 *  specifically, transpose the multiplier and use a modified
 *  multiplication algorithm with the transposed multiplier.
 *
 *  Set *err to EINVAL if nRows or nCols <= 0, to ENOMEM if not enough
 *  memory.
 */
SmartMulMatrix *newSmartMulMatrix(int nRows, int nCols, int *err)
{
  SmartMulMatrixImpl *matrix = (SmartMulMatrixImpl *)newDenseMatrix(nRows, nCols, err);
  if (*err == EINVAL || *err == EDOM) return NULL;
  
  matrix->fns = (MatrixFns *)getSmartMulMatrixFns();
  return (SmartMulMatrix *)matrix;
}

static void patchSmartMulMatrixFns(void)
{
  if (!isInit) {
    const DenseMatrixFns *fns = getDenseMatrixFns();
    smartMulMatrixFns.free = fns->free;
    smartMulMatrixFns.getNRows = fns->getNRows;
    smartMulMatrixFns.getNCols = fns->getNCols;
    smartMulMatrixFns.getElement = fns->getElement;
    smartMulMatrixFns.setElement = fns->setElement;
    smartMulMatrixFns.transpose = fns->transpose;
    isInit = true;
  }
}

/** Return implementation of functions for a smart multiplication
 *  matrix; these functions can be used by sub-classes to inherit
 *  behavior from this class.
 */
const SmartMulMatrixFns *
getSmartMulMatrixFns(void)
{
  patchSmartMulMatrixFns();
  return &smartMulMatrixFns;
}
