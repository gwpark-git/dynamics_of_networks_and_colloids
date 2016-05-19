
#ifndef REPULSIVE_BROWNIAN_H
#define REPULSIVE_BROWNIAN_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"

MKL_LONG main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition);

#endif
