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

    double calculateValue(const Walker *walker);

};

}
