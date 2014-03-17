#pragma once

#include <QMC2.h>
#include "HF/src/libs/system/system.h"

namespace QMC2
{

class HFOrbitals : public Orbitals
{
public:
    HFOrbitals(hf::System *system, mat corePositions);

    // Orbitals interface
public:
    void set_qnum_indie_terms(Walker *walker, int i) {}
    double phi(const Walker *walker, int particle, int q_num);
    double del_phi(const Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const Walker *walker, int particle, int q_num);

    hf::System* m_system;

    mat m_corePositions;
};

}
