#ifndef EXPANDEDBASISFUNCTIONS_H
#define EXPANDEDBASISFUNCTIONS_H

#include "BasisFunctions.h"

#include <armadillo>
#include <vector>

class expandedBasisFunctions : public BasisFunctions
{
public:
    expandedBasisFunctions(std::vector<BasisFunctions*> BF, arma::vec coeffs) : BF(BF), coeffs(coeffs), n(coeffs.n_elem) {}
private:
    std::vector<BasisFunctions*> BF;
    arma::vec coeffs;
    int n;

public:
    double eval(const Walker *walker, int i) {
        int k = 0;
        double val = 0;

        while (k != n) {
            val += coeffs(k)*BF.at(k)->eval(walker, i);
            k++;
        }

        return val;
    }
};

#endif // EXPANDEDBASISFUNCTIONS_H
