
#ifndef STOCHASTIC_HEUR_FLOWERS_H
#define STOCHASTIC_HEUR_FLOWERS_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/association.h"
#include "../lib/handle_association.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"


MKL_LONG stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition);


#endif
