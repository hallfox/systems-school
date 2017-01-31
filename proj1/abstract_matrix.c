#include "abstract_matrix.h"

#include <errno.h>
#include <stdlib.h>

static const char *getKlass(const Matrix *this, int *err)
{
  return "abstractMatrix";
}

static void freeAbstractMatrix(Matrix *this, int *err)
{
  free(this);
}

static void transpose(const Matrix *this, Matrix *result, int *err)
{
  // Check dimensions: MxN -> NxM
  const int this_m = this->fns->getNRows(this, err);
  if (*err == EINVAL) return;
  const int this_n = this->fns->getNCols(this, err);
  if (*err == EINVAL) return;
  const int result_n = this->fns->getNRows(this, err);
  if (*err == EINVAL) return;
  const int result_m = this->fns->getNCols(this, err);
  if (*err == EINVAL) return;
  if (!(this_m == result_m && this_n == result_n)) {
    *err = EDOM;
    return;
  }

  // Transpose: Res[r][c] = This[c][r]
  for (int r = 0; r < result_n; r++) {
    for (int c = 0; c < result_m; c++) {
      MatrixBaseType x = this->fns->getElement(this, c, r, err);
      if (*err == EDOM || *err == EINVAL) return;
      result->fns->setElement(result, r, c, x, err);
      if (*err == EDOM || *err == EINVAL) return;
    }
  }
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

  // Do the multiplication
  for (int pr_r = 0; pr_r < pr_m; pr_r++) {
    for (int pr_c = 0; pr_c < pr_p; pr_c++) {
      // Pr[r][c] <- Sum_i This[r][i]*That[i][c]
      MatrixBaseType res = 0;
      for (int i = 0; i < this_n; i++) {
	MatrixBaseType a = this->fns->getElement(this, pr_r, i, err);
	if (*err == EINVAL || *err == EDOM) return;
	MatrixBaseType b = multiplier->fns->getElement(this, i, pr_c, err);
	if (*err == EINVAL || *err == EDOM) return;
	res += a*b;
      }
      product->fns->setElement(product, pr_r, pr_c, res, err);
      if (*err == EINVAL || *err == EDOM) return;
    }
  }
}

static MatrixFns abstractMatrixFns = {
  .getKlass = getKlass,
  .free = freeAbstractMatrix,
  .transpose = transpose,
  .mul = mul,
};

/** Return implementation of functions for an abstract matrix; these are
 *  functions which can be implemented using only other matrix functions,
 *  independent of the actual implementation of the matrix.
 */
const MatrixFns *
getAbstractMatrixFns(void)
{
  return &abstractMatrixFns;
}
