#pragma once

/*	Reference: Jing, RJ., Moreno-Maza, M., Talaashrafi, D. (2020). 
	Complexity Estimates for Fourier-Motzkin Elimination. 
	In: Boulier, F., England, M., Sadykov, T.M., Vorozhtsov, E.V. (eds) 
	Computer Algebra in Scientific Computing. CASC 2020. 
	Lecture Notes in Computer Science(), vol 12291. Springer, Cham. 
	https://doi.org/10.1007/978-3-030-60026-6_16
*/

#include "matrixOperations.h"

void partition(std::vector<std::vector<mpq_class>> &J,
				const std::vector<mpq_class> &A, 
				std::vector<std::vector<mpq_class>> &Jplus, 
				std::vector<std::vector<mpq_class>> &Jzero,
				std::vector<std::vector<mpq_class>> &Jminus);

// adjacency test based on https://arxiv.org/pdf/1811.01510 Alg. 1.
bool adjacencyTest(std::vector<std::vector<mpq_class>> &Ak, 
	std::vector<mpq_class> &Jplusp, std::vector<mpq_class> &Jminusn);

std::vector<std::vector<mpq_class>> ddMethod(std::vector<std::vector<mpq_class>> &A);



