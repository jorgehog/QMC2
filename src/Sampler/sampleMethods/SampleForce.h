#ifndef SAMPLEFORCE_H
#define SAMPLEFORCE_H

#include "../Sampler.h"

class SampleForce : public Sampler
{
public:

    SampleForce(double *R, int n_p);

private:

    double * R;
    int n_p;

    double calculateValue(const Walker *walker);

};

#endif // SAMPLEFORCE_H
