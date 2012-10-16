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

    virtual void initialize(Walker* walker);
    virtual void get_spatial_grad(Walker* walker, int particle) const;
    virtual void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int i);
    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const;

    virtual double get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const;
    virtual double get_spatial_lapl_sum(const Walker* walker) const;
    virtual double get_spatial_wf(const Walker* walker);

    virtual void copy_walker(const Walker* parent, Walker* child) const;
    virtual void reset_walker_ISCF(const Walker* walker_pre, Walker* walker_post, int particle) const;

};

#endif	/* FERMIONS_H */

