#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "obs_stat/Constants.h"

//For the information gain
#define MIN_LOG_VAL 1e-3
#define MIN_POSITIVE_VALUE 1e-8
#define MAX_NOBS_FOR_IG 1e5
#define OCC_THRESHOLD 0.8

//For the costmap 
#define FREE 0xFF
#define UNKNOWN 0x80
#define OCCUPIED 0x00
#define CMAP_LETHAL 100
#define CMAP_UNKNOWN -1
#define CMAP_MAXVALID 98

#endif