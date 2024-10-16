#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cassert>
#include <gmp.h>
#include <gmpxx.h>

/* a matrix is implemented as a vector of vectors of mpq_t elements. */

/* matrix A has m rows and n columns. */
void makeRandomMatrix(std::vector<std::vector<mpq_class>> &A, int m, int n);

/* print a matrix */
void printMatrix(const std::vector<std::vector<mpq_class>> &A);

/* convert a matrix to row-reduced echelon form
   based on: https://rosettacode.org/wiki/Reduced_row_echelon_form#C
   @param order keeps track of the row swap.
*/
void matrixToRREF(std::vector<std::vector<mpq_class>> &A, std::vector<int> &order);

/* returns the index of first row that all elements are zero. this means that
   all rows before this index are independent.
*/
int findFirstZeroRow(std::vector<std::vector<mpq_class>> &RA);

/* returns the rank of a matrix 
   based on https://cp-algorithms.com/linear_algebra/rank-matrix.html
*/
int computeRank(std::vector<std::vector<mpq_class>> A);

/* computing matrix inverse via LU factorization.
   based on: https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
   P.size() = A.size() + 1
 */
void LUPDecompose(std::vector<std::vector<mpq_class>> &A, std::vector<int> &P);

void LUPInvert(const std::vector<std::vector<mpq_class>>& A, const std::vector<int>& P, 
				std::vector<std::vector<mpq_class>>& IA);

void rowColumnMultiply(mpq_class &prod, const std::vector<mpq_class> &row, 
						const std::vector<mpq_class> &col);

// prod = mat partition* col. 
void matrixVectorMultiply(std::vector<mpq_class> &prod, 
	std::vector<std::vector<mpq_class>> &mat,
	const std::vector<mpq_class> &col);

/* adjust the output matrix to have only 0 (vertex) or 1 (ray) at the last row. */
void lastRowZeroAndOne(std::vector<std::vector<mpq_class>> &R0);

