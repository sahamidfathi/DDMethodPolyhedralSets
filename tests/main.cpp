#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include "../include/matrixOperations.h"
#include "../include/ddMethod.h"

// returns true if VR1 should be placed BEFORE VR2 (in REVLEX ORDER)
// returns false if VR1 should be placed AFTER VR2 (in REVLEX ORDER)
bool VRCompare(std::vector<mpq_class> VR1, std::vector<mpq_class> VR2) {
    for (int i = 0; i < VR1.size(); i++) {
        // we use equal first because it is way faster than compare
        int equalResult = mpq_equal(VR1.at(i).get_mpq_t(), VR2.at(i).get_mpq_t());
        // equalResult == 0 means that the two mpqs are NOT EQUAL
        if (equalResult == 0) {
			int cmpResult = mpq_cmp(VR1.at(i).get_mpq_t(), VR2.at(i).get_mpq_t());
			// cmpResult is negative -> input1 < input2
			if (cmpResult < 0) {
				return false;
			}
			else {
				return true;
			}
		}
        // if they ARE equal we simply continue looping
	}
    // if we have reached the end and all mpqs were found to be equal, it does not matter what we return
    return true;
}

// sorting vertices and rays.
void sortVR(std::vector<std::vector<mpq_class>> &VR) {
    std::sort(VR.begin(), VR.end(), VRCompare);
}

// transpose matrix (only for comparing results).
std::vector<std::vector<mpq_class>> transpose(std::vector<std::vector<mpq_class>> A) {
	std::vector<std::vector<mpq_class>> At(A[0].size(), 
											std::vector<mpq_class>(A.size()));
	for (int i = 0; i < A.size(); ++i)	{
		for (int j = 0; j < A[0].size(); ++j) {
			At[i][j] = A[j][i];
		}
	}
	return At;
}

std::vector<std::vector<mpq_class>> readMatrixFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<mpq_class>> matrix;
    
    if (!file.is_open()) {
        std::cerr << "Could not open the file " << filename << std::endl;
        return matrix;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<mpq_class> row;
        mpq_class value;
        while (iss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

void writeMatrixToFile(const std::string& filename, const std::vector<std::vector<mpq_class>>& matrix) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Could not open the file " << filename << " for writing." << std::endl;
		return;
	}
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            file << value.get_str() << " ";
        }
        file << std::endl;
    }
    file.close();
}

int main(int argc, char **argv) {
	
	if (argc != 3) {
		std::cerr << "Usage: " << 
		argv[0] << " <input_path> <output_path>" << std::endl;
		return 1;
	}

    std::string filenameIn = argv[1];
    std::string filenameOut = argv[2];
	
	std::vector<std::vector<mpq_class>> A = readMatrixFromFile(filenameIn);
	
	std::vector<std::vector<mpq_class>> resultVecs = ddMethod(A);
	sortVR(resultVecs);
	
	writeMatrixToFile(filenameOut, resultVecs);

	return 0;
}



