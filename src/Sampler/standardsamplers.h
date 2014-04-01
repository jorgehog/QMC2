#pragma once

#include "Sampler.h"

#include "../Walker/Walker.h"

namespace QMC2
{

class SampleRadius : public Sampler
{
    // Sampler interface
public:
    void push_values(Walker *walker)
    {
        for (const double & r : walker->abs_r)
        {
            push_value(r);
        }
    }
};

class CalculateAndSampleRadius : public Sampler
{
    // Sampler interface
public:
    void push_values(Walker *walker)
    {
        walker->calc_r_i();

        for (const double & r : walker->abs_r)
        {
            push_value(r);
        }
    }
};

class SampleRadiusSquared : public Sampler
{
    // Sampler interface
public:
    void push_values(Walker *walker)
    {
        for (const double & r : walker->r2)
        {
            push_value(r);
        }
    }
};

class CalculateAndSampleRadiusSquared : public Sampler
{
    // Sampler interface
public:
    void push_values(Walker *walker)
    {
        walker->calc_r_i2();

        for (const double & r : walker->r2)
        {
            push_value(r);
        }
    }
};

class SampleRelativeDistance : public Sampler
{
    // Sampler interface
public:
    void push_values(Walker *walker)
    {
        for (uint i = 0; i < walker->r_rel.n_cols; ++i)
        {
            for (uint j = i + 1; j < walker->r_rel.n_cols; ++j)
            {
                push_value(walker->r_rel(i, j));
            }
        }
    }
};

}
