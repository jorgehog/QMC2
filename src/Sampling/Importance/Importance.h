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

    virtual void get_necessities(Walker* walker);
    virtual void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle);
    virtual void calculate_energy_necessities_CF(Walker* walker) const;

    virtual double get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const;

    virtual void copy_walker(const Walker* parent, Walker* child) const;
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

};

#endif	/* IMPORTANCE_H */

