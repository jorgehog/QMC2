/* 
 * File:   Fermions.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:44 PM
 */

#ifndef FERMIONS_H
#define	FERMIONS_H

class Fermions : public System {
protected:

    arma::rowvec I;

    void make_merged_inv(Walker* walker);
    void update_inverse(const Walker* walker_old, Walker* walker_new, int particle);


public:
    Fermions(GeneralParams &, Orbitals*);

    virtual void get_spatial_grad(Walker* walker, int particle) const;
    virtual void get_spatial_grad_full(Walker* walker) const;
    virtual double get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const;
    virtual double get_spatial_lapl_sum(const Walker* walker) const;

    virtual void copy_walker(const Walker* parent, Walker* child) const {
        child->inv = parent->inv;
    }

    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
        using namespace arma;
        walker_pre->inv(span(0, n2 - 1), span(start, end)) = walker_post->inv(span(0, n2 - 1), span(start, end));
    };

    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
        using namespace arma;
        walker_post->inv(span(0, n2 - 1), span(start, end)) = walker_pre->inv(span(0, n2 - 1), span(start, end));
    }

    virtual double get_spatial_wf(const Walker* walker) {
        using namespace arma;
        return det(walker->phi(span(0, n2 - 1), span())) * det(walker->phi(span(n2, n_p - 1), span()));
    }

    virtual void initialize(Walker* walker) {
        make_merged_inv(walker);
    }

    virtual void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int i) {
        update_inverse(walker_old, walker_new, i);
    }

};

#endif	/* FERMIONS_H */

