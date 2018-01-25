IMPORT $ AS LR;
IMPORT LR.Types AS LR_Types;
IMPORT Std.BLAS.Types AS Types;

/**
 * Create a (row or column) vector from a general matrix using a function
 * such as max, min, average, max absolute value, et cetera.
 * <table>
 *  <tr>
 *    <th>Operation</th>
 *    <th>Description</th>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Max_Val</td>
 *    <td>Maximum value in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Min_Val</td>
 *    <td>MInimum value in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Ave_Val</td>
 *    <td>Average of the values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Sum_Val</td>
 *    <td>Sum of the values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Ave_Absl</td>
 *    <td>Average of the absolute value of the values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Sum_Abs</td>
 *    <td>Sum of the absolute value of the values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Sum_Sqs</td>
 *    <td>Sum of the value squared of the values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Zeros</td>
 *    <td>Number of zero values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Non_Zeros</td>
 *    <td>Number of non-zero values in the row or column.</td>
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Max_Pos</td>
 *    <td>Position of the maximum value, zero based.
 *  </tr>
 *  <tr>
 *    <td>Matrix_Op.Min_Pos</td>
 *    <td>Position of the minimum value, zero based.
 *  </tr>
 * </table>
 * @param M number of rows
 * @param N number of columns
 * @param A the MxN matrix, column major form
 * @param row_vector flag to produce a row vector of the column values if
 * true, or a column vector of the row values if false
 * @param op the operation enumeration value
 * @returns the row or column vector, 1xN or Mx1
 */
EXPORT Types.matrix_t gemagg(Types.dimension_t M, Types.dimension_t N,
                             Types.matrix_t A, BOOLEAN row_vector,
                             LR_Types.Matrix_Op op) := BEGINC++
  const uint8_t op_max_val = 1;
  const uint8_t op_min_val = 2;
  const uint8_t op_ave_val = 3;
  const uint8_t op_sum_val = 4;
  const uint8_t op_ave_abs = 5;
  const uint8_t op_sum_abs = 6;
  const uint8_t op_sum_sqs = 7;
  const uint8_t op_zeros = 8;
  const uint8_t op_nonzeros = 9;
  const uint8_t op_max_pos = 10;
  const uint8_t op_min_pos = 11;
  if (lenA < m*n*sizeof(double)) rtlFail(0, "gemo: MxN larger than matrix");
  size32_t step = (row_vector) ? 1 : m;
  size32_t stop = (row_vector) ? m : m*n;
  size32_t pos_step = (row_vector) ? m : 1;
  size32_t answers = (row_vector) ? n : m;
  __lenResult = answers * sizeof(double);
  __result = rtlMalloc(__lenResult);
  __isAllResult = false;
  double* rslt = (double*) __result;
  double* mat = (double*) a;
  size32_t pos = 0;
  for (size32_t i=0; i<answers; i++) {
    switch(op) {
      case op_max_val:
      case op_min_val:
        rslt[i] = mat[pos + (i*pos_step)];
        break;
      case op_max_pos:
      case op_min_pos:
        rslt[i] = pos + (i*pos_step);
        break;
      default:
        rslt[i] = 0;
        break;
    }
  }
  for (size32_t i=0; i<answers; i++) {
    switch (op) {
      case op_max_val:
        for (size32_t j=0; j<stop; j+=step) if (rslt[i]<mat[pos+j]) rslt[i]=mat[pos+j];
        break;
      case op_min_val:
        for (size32_t j=0; j<stop; j+=step) if (rslt[i]>mat[pos+j]) rslt[i]=mat[pos+j];
        break;
      case op_ave_val:
        for (size32_t j=0; j<stop; j+=step) rslt[i] += mat[pos+j];
        rslt[i] = rslt[i] / ((row_vector) ? m : n);
        break;
      case op_sum_val:
        for (size32_t j=0; j<stop; j+=step) rslt[i] += mat[pos+j];
        break;
      case op_ave_abs:
        for (size32_t j=0; j<stop; j+=step) rslt[i] += abs(mat[pos+j]);
        rslt[i] = rslt[i] / ((row_vector) ? m : n);
        break;
      case op_sum_abs:
        for (size32_t j=0; j<stop; j+=step) rslt[i] += abs(mat[pos+j]);
        break;
      case op_sum_sqs:
        for (size32_t j=0; j<stop; j+=step) rslt[i] += mat[pos+j]*mat[pos+j];
        break;
      case op_zeros:
        for (size32_t j=0; j<stop; j+=step) if(mat[pos+j]==0) rslt[i]++;
        break;
      case op_nonzeros:
        for (size32_t j=0; j<stop; j+=step) if(mat[pos+j]!=0) rslt[i]++;
        break;
      case op_max_pos:
        for (size32_t j=0; j<stop; j+=step) if(mat[(size32_t)rslt[i]]<mat[pos+j]) rslt[i]=pos+j;
        break;
      case op_min_pos:
        for (size32_t j=0; j<stop; j+=step) if(mat[(size32_t)rslt[i]]>mat[pos+j]) rslt[i]=pos+j;
        break;
    }
    pos += pos_step;
  }
ENDC++;
