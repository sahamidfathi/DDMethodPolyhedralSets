#include "../include/matrixOperations.h"

void makeRandomMatrix(std::vector<std::vector<mpq_class>> &A, int m, int n) {
	A.resize(m);
	for (int i = 0; i < m; ++i) {
		A[i].resize(n);
		for (int j = 0; j < n; ++j) {
			A[i][j] = mpq_class(i + 1, j + 2);
		}
	}
}

void printMatrix(const std::vector<std::vector<mpq_class>> &A) {
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
			//mpq_out_str(stdout, 10, A[i][j].get_mpq_t ());
			mpq_out_str(stdout, 10, A[i][j].get_mpq_t());
			std::cout << ",\t";
		}
		std::cout << '\n';
	}
}

void matrixToRREF(std::vector<std::vector<mpq_class>> &A, std::vector<int> &order) {
	int lead;
	int rix, iix;
	mpq_class leadVal;
	int rowCount = A.size();
	int colCount = A[0].size();

	lead = 0;
	for (rix = 0; rix < rowCount; ++rix) {
		if (lead >= colCount)
			return;

		iix = rix;
		while (0 == A[iix][lead]) {
			iix++;
			if (iix == rowCount) {
				iix = rix;
				lead++;
				if (lead == colCount)
					return;
			}
		}

		if (iix != rix) {
			mpq_class tmp;
			for (int ix = 0; ix < colCount; ++ix) {
				tmp = A[iix][ix];
				A[iix][ix] = A[rix][ix];
				A[rix][ix] = tmp;
			}
			// recording order change
			int tmpindex = order[iix];
			order[iix] = order[rix];
			order[rix] = tmpindex;
		}

		leadVal = A[rix][lead];
		if (leadVal) {
			for (int ix = 0; ix < colCount; ++ix)
				A[rix][ix] /= leadVal;
		}

		for (iix = 0; iix < rowCount; ++iix) {
			if (iix != rix) {
				leadVal = A[iix][lead];
				for (int ix = 0; ix < colCount; ++ix)
				A[iix][ix] += (-leadVal * A[rix][ix]);
			}
		}
		lead++;
	}
}

int findFirstZeroRow(std::vector<std::vector<mpq_class>> &RA) {
	bool allZero;
	for (int i = 0; i < RA.size(); ++i) {
		allZero = true;
		for (int j = 0; j < RA[1].size(); ++j) {
			if (RA[i][j]) {
				allZero = false;
				break;
			}
		}
		if (allZero) {
			return i;
		}
	}
	return RA.size(); // all rows are linearly independent.
}

int computeRank(std::vector<std::vector<mpq_class>> A) {
	int m = A.size();
	int n = A[0].size();

	int rank = 0;
	std::vector<bool> rowSelected(m, false);
	for (int i = 0; i < n; ++i) {
		int j;
		for (j = 0; j < m; ++j) {
			if (!rowSelected[j] && A[j][i])
			break;
		}

		if (j != m) {
			++rank;
			rowSelected[j] = true;
			for (int p = i + 1; p < n; ++p) {
				A[j][p] /= A[j][i];
			}
			for (int k = 0; k < m; ++k) {
				if (k != j && A[k][i]) {
					for (int p = i + 1; p < n; ++p) {
						A[k][p] -= A[j][p] * A[k][i];
					}
				}
			}
		}
	}
	return rank;
}

void LUPDecompose(std::vector<std::vector<mpq_class>> &A, std::vector<int> &P) {

	assert(A.size() == A[0].size() && "Matrix must be square.");
	
	int dim = A.size(); // number of rows or columns.
	
	int i, j, k, imax; 
	mpq_class maxA, absA;

	for (i = 0; i <= dim; ++i)
		P[i] = i; //unit permutation matrix, P[dim] initialized with dim

	for (i = 0; i < dim; ++i) {
		maxA = 0;
		imax = i;

		for (k = i; k < dim; ++k)
			if ((absA = abs(A[k][i])) > maxA) { 
				maxA = absA;
				imax = k;
			}

		if (imax != i) {
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			std::swap(A[i], A[imax]);

			//counting pivots starting from dim (for determinant)
			P[dim]++;
		}

		for (j = i + 1; j < dim; ++j) {
			A[j][i] /= A[i][i];

			for (k = i + 1; k < dim; ++k)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}
}

void LUPInvert(const std::vector<std::vector<mpq_class>>& A, const std::vector<int>& P, 
				std::vector<std::vector<mpq_class>>& IA) {
	int dim = A.size();
	for (int j = 0; j < dim; ++j) {
		for (int i = 0; i < dim; ++i) {
			IA[i][j] = (P[i] == j) ? 1 : 0;
			for (int k = 0; k < i; ++k)
				IA[i][j] -= A[i][k] * IA[k][j];
		}

		for (int i = dim - 1; i >= 0; --i) {
			for (int k = i + 1; k < dim; ++k)
				IA[i][j] -= A[i][k] * IA[k][j];

			IA[i][j] /= A[i][i];
		}
	}
}

void rowColumnMultiply(mpq_class &prod, const std::vector<mpq_class> &row, 
						const std::vector<mpq_class> &col) {
	prod = 0;
	for (int i = 0; i < row.size(); ++i) {
		prod = prod + row[i] * col[i];
	}
}

void matrixVectorMultiply(std::vector<mpq_class> &prod, 
	std::vector<std::vector<mpq_class>> &mat,
	const std::vector<mpq_class> &col) {
	
	mpq_class tmp;
	for (int i = 0; i < mat.size(); ++i) {
		tmp = 0;
		for (int j = 0; j < mat[0].size(); ++j) {
			tmp = tmp + mat[i][j] * col[j];
		}
		prod[i] = tmp;
	}	
}

void lastRowZeroAndOne(std::vector<std::vector<mpq_class>> &R0) {

	if (R0.size() < 1)
		return;

	for (int j = 0; j < R0[0].size(); ++j) {
	
		if (R0[R0.size() - 1][j] == 0)
			continue;
			
		mpz_class denom = R0[R0.size() - 1][j].get_den();
		mpz_class num = R0[R0.size() - 1][j].get_num();
		for (int i = 0; i < R0.size(); ++i) {
			R0[i][j] *= denom;
			R0[i][j] /= num;	
		}
	}
}

