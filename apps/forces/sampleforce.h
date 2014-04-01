#pragma once


#include <QMC2.h>

namespace QMC2
{

class SampleForce : public Sampler
{
public:

    SampleForce(double *R, int n_p);

private:

    double * R;
    int n_p;



    // Sampler interface
public:
    void push_values(const Walker *walker);
};

}
