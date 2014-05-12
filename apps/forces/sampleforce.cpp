#include "sampleforce.h"


using namespace QMC2;

SampleForce::SampleForce(double * R, int n_p) :
    Sampler(), n_p(n_p)
{
    this->R = R;
}

void SampleForce::push_values(const Walker *walker)
{

    double Rhalf = *R/2;

    double force = 0;

    double com_corr, shared;
    double quarterR2 = Rhalf*Rhalf;

    for (int i = 0; i < n_p; i++) {
/*
        EXPLAIRNATION

        |R/2 +- r| = sqrt|(R/2 - rx)^2 + ry^2 + rz^2|
                   = sqrt|(R/2)^2 + |r|^2 +-R*rx|

            shared == (R/2)^2 + |r|^2
          com_corr == R*rx
        |R/2 +- r| = sqrt|shared +- com_corr| == rp(lus) or rm(inus)
*/

        double rx = walker->r(i, 0);

        shared = walker->get_r_i2(i)+ quarterR2;
        com_corr = (*R)*rx;

        double rp = sqrt(shared + com_corr);
        double rm = sqrt(shared - com_corr);

        double thresh = 0.01;
        if ((rp*rp < thresh*thresh) || (rm*rm < thresh*thresh)){
            //std::cout << "WARNING IN FORCE CALC: " << rp << "  " << rm << std::endl;
            return;
        }

        force += (Rhalf + rx)/(rp*rp*rp) + (Rhalf - rx)/(rm*rm*rm);

    }

    force -= 1. /(*R*(*R));

    push_value(0.5*force);
}
