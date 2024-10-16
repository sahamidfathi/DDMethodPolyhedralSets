#include "../include/ddMethod.h"

void partition(std::vector<std::vector<mpq_class>> &J,
				const std::vector<mpq_class> &A, 
				std::vector<std::vector<mpq_class>> &Jplus, 
				std::vector<std::vector<mpq_class>> &Jzero,
				std::vector<std::vector<mpq_class>> &Jminus) {
	
	mpq_class prod;
	
	for (int i = 0; i < J.size(); ++i) {
		rowColumnMultiply(prod, A, J[i]);
		if (prod > 0) {
			Jplus.emplace_back(J[i]);
		} else if (prod == 0) {
			Jzero.emplace_back(J[i]);
		} else {
			Jminus.emplace_back(J[i]);
		}
	}
}

bool adjacencyTest(std::vector<std::vector<mpq_class>> &Ak, 
	std::vector<mpq_class> &Jplusp, std::vector<mpq_class> &Jminusn) {
	
	std::vector<mpq_class> s0(Ak.size());
	std::vector<mpq_class> s1(Ak.size());	
	
	matrixVectorMultiply(s0, Ak, Jplusp);
	matrixVectorMultiply(s1, Ak, Jminusn);

	std::vector<std::vector<mpq_class>> subMatAk;
	for (int i = 0; i < s0.size(); ++i) {
		if ((s0[i] == 0) && (s1[i] == 0))
			subMatAk.emplace_back(Ak[i]);
	}
	
	// find rank of subMatAk
	int rank = 0;
	
	if (subMatAk.size() > 0) 
		rank = computeRank(subMatAk);
	
	return (rank == (Ak[0].size() - 2));
}

std::vector<std::vector<mpq_class>> ddMethod(std::vector<std::vector<mpq_class>> &A) {

	printMatrix(A);
	
	// order holds the order of the rows of the matrix
	std::vector<int> order(A.size());
	for (int i = 0; i < A.size(); ++i) 
		order[i] = i;
	
//	std::cout << "order before RREF is: " << std::endl;
//	for (int i = 0; i < A.size(); ++i) 
//		std::cout << order[i] << std::endl;
	
	std::vector<std::vector<mpq_class>> RA(A);
	matrixToRREF(RA, order);
//	std::cout << "RREF is: " << std::endl;
//	printMatrix(RA);
	
//	std::cout << "order after RREF is: " << std::endl;
//	for (int i = 0; i < A.size(); ++i) 
//		std::cout << order[i] << std::endl;

	// index of first row with all zeros (start of independent rows).
	int zrow = findFirstZeroRow(RA);
//	std::cout << "first row with all zero in RREF is: " << zrow << std::endl;
		
	// Ak is a sub-matrix of independent rows of A.
	std::vector<std::vector<mpq_class>> Ak(zrow, std::vector<mpq_class>(A[0].size()));
	for (int i = 0; i < zrow; ++i) {
		Ak[i] = A[order[i]];
	}
//	std::cout << "Ak is: " << std::endl;
//	printMatrix(Ak);
	
	// vector of indices of non-independent rows of A 
	std::vector<int> depIdx(A.size() - zrow);
	
	for (int i = zrow; i < A.size(); ++i) {
		depIdx[i - zrow] = order[i];
	}
//	std::cout << "indices of non-ind. rows: "<< std::endl;
//	for (int i = 0; i < depIdx.size(); ++i)
//		std::cout << depIdx[i] << std::endl;
	
	std::vector<std::vector<mpq_class>> origAk = Ak;
	
	// R' is inverse of Ak.
	std::vector<int> P(Ak.size() + 1);
	LUPDecompose(Ak, P);
	
	std::vector<std::vector<mpq_class>> IAk(Ak.size(), 
										std::vector<mpq_class>(Ak[0].size()));
	LUPInvert(Ak, P, IAk);
	
//	std::cout << "Ak^(-1) is: " << std::endl;
//	printMatrix(IAk);
	
	// J is a vector of columns of R' or IAk (IAk transpose)
	std::vector<std::vector<mpq_class>> J(IAk[0].size(), 
											std::vector<mpq_class>(Ak.size()));
	for (int i = 0; i < J.size(); ++i) {
		for (int j = 0; j < J[0].size(); ++j) {
			J[i][j] = IAk[j][i];
		}
	}
	std::vector<std::vector<mpq_class>> R0End = IAk;
	
	// loop over indices of dependent rows of A.
	for (int i = 0; i < depIdx.size(); ++i) {
		// select A[depIdx[i]]
		
		// set R' (= R0) to be an empty matrix
		std::vector<std::vector<mpq_class>> R0(origAk[0].size()); // = IAk;
	
		// partition 
		std::vector<std::vector<mpq_class>> Jplus;
		std::vector<std::vector<mpq_class>> Jzero;
		std::vector<std::vector<mpq_class>> Jminus;

		partition(J, A[depIdx[i]], Jplus, Jzero, Jminus);

		// add vectors in Jplus and Jzero as columns to R0.
		// (alternatively, we can add vectors of Jplus and Jzero as rows to R0 and 
		// transpose R0 at the end.)
		for (int ii = 0; ii < Jplus.size(); ++ii) {
			for (int j = 0; j < Jplus[ii].size(); ++j) {
				R0[j].emplace_back(Jplus[ii][j]);
			}
		}
		for (int ii = 0; ii < Jzero.size(); ++ii) {
			for (int j = 0; j < Jzero[ii].size(); ++j) {
				R0[j].emplace_back(Jzero[ii][j]);
			}
		}
		
//		std::cout << "R0 after is " << std::endl;
//		printMatrix(R0);
		
		// do adjacency test.
		for (int p = 0; p < Jplus.size(); ++p) {
			for (int n = 0; n < Jminus.size(); ++n) {
				if (adjacencyTest(origAk, Jplus[p], Jminus[n])) { 
					// dotProduct(A[depIdx[i]], Jplus[p]);
					// rNew = (Ai*rp)*rn - (Ai*rn)*rp
					
					mpq_class Ai_rp;
					mpq_class Ai_rn;
					
					rowColumnMultiply(Ai_rp, A[depIdx[i]], Jplus[p]);
					rowColumnMultiply(Ai_rn, A[depIdx[i]],  Jminus[n]);					
					
					for (int k = 0; k < A[0].size(); ++k) {
						R0[k].emplace_back(Ai_rp * Jminus[n][k] - Ai_rn * Jplus[p][k]);
					}
				}
			}
		}
		
		std::vector<std::vector<mpq_class>> Jnew(R0[0].size(), 
					std::vector<mpq_class>(Ak.size()));
		for (int i = 0; i < R0[0].size(); ++i) {
			for (int j = 0; j < R0.size(); ++j) {
				Jnew[i][j] = R0[j][i];
			}
		}
		
		J = Jnew;
		
		origAk.emplace_back(A[depIdx[i]]);
		
		R0End = R0;
	}
	
	lastRowZeroAndOne(R0End);
	
    // transpose R0.
    std::vector<std::vector<mpq_class>> R0T(R0End[0].size(), 
                                            std::vector<mpq_class>(R0End.size()));
    for (int i = 0; i < R0T.size(); ++i) {
        for (int j = 0; j < R0T[0].size(); ++j) {
            R0T[i][j] = R0End[j][i];
        }
    }
    
    R0End = R0T;

//    std::cout << "R0End is " << std::endl;
//    printMatrix(R0End);

	return R0End;
}







