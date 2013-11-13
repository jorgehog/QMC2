#ifndef OXYGEN3_21G_H
#define OXYGEN3_21G_H

#include "../gaussianfitted.h"

struct GeneralParams;
struct VariationalParams;

#include "../../../Walker/Walker.h"
#include "../../../BasisFunctions/BasisFunctions.h"

class Oxygen3_21G : public GaussianFitted
{
public:
    Oxygen3_21G();

//    double phi(const Walker *walker, int particle, int q_num) {

//        for (double* ef : expFactors) {
//            std::cout << *ef << std::endl;
//        }

//        return basis_functions[q_num]->eval(walker, particle);
//    }

    double del_phi(const Walker *walker, int particle, int q_num, int d) {
        return num_diff(walker, particle, q_num, d);
    }

    double lapl_phi(const Walker *walker, int particle, int q_num){
        return num_ddiff(walker, particle, q_num);
    }
};

#endif // OXYGEN3_21G_H
