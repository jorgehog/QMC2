/* 
 * File:   Importance.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:43 PM
 */

#ifndef IMPORTANCE_H
#define	IMPORTANCE_H

class Importance : public Sampling {
public:

    Importance(GeneralParams &);

    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const;

    virtual void get_necessities(Walker* walker) {
        qmc->get_gradients(walker);
        qmc->get_QF(walker);
    }

    virtual void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
        qmc->get_gradients(walker_pre, walker_post, particle);
        qmc->get_QF(walker_post);
    }

    virtual void calculate_energy_necessities(Walker* walker) const {

    }

    virtual void copy_walker(const Walker* parent, Walker* child) const {
        child->jast_grad = parent->jast_grad;
        child->spatial_grad = parent->spatial_grad;
        child->qforce = parent->qforce;
    }
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

};

#endif	/* IMPORTANCE_H */

