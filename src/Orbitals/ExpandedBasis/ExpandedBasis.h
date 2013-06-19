/* 
 * File:   ExpandedBasis.h
 * Author: jorgmeister
 *
 * Created on October 9, 2012, 5:09 PM
 */

#ifndef EXPANDEDBASIS_H
#define	EXPANDEDBASIS_H

#include "../Orbitals.h"
struct GeneralParams;

class ExpandedBasis : public Orbitals {
public:
    ExpandedBasis(GeneralParams & gp, Orbitals* basis, int m, std::string coeffPath);

    double phi(const Walker* walker, int particle, int q_num);
    double del_phi(const Walker* walker, int particle, int q_num, int d);
    double lapl_phi(const Walker* walker, int particle, int q_num);

    void set_qnum_indie_terms(Walker* walker, int i) {
        basis->set_qnum_indie_terms(walker, i);
    }
    
protected:
    int basis_size;

    arma::mat coeffs;

    Orbitals* basis;

};

#endif	/* EXPANDEDBASIS_H */

