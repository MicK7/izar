#ifndef IZARDEFINES_H
#define IZARDEFINES_H

#include "vtkSetGet.h"

/* Naming conventions */

// Gas model
#define GAMMA "Gamma"
#define RGAS "Rgas"
#define CP "Cp"

// Gas mixing names
#define FAR "FAR"
#define WAR "WAR"
#define RFUEL "Rfuel"
#define RWATER "Rwater"
#define RMIX "Rmix"

// Cp(T) polynomial model names
#define CPAIR "CpPolyAir"
#define CPFUEL "CpPolyFuel"
#define CPWATER "CpPolyWater"
#define CPMIX "CpPolyMix"
#define TS_LOWER_BOUND 200.0
#define TS_UPPER_BOUND 2000.0
#define TS_PRESSURE_EFFECT_LIMIT  1800.0

// Turbomachinery stuff
#define OMEGA "Omega"
#define ZSECTOR "Zsector"
#define ROW_NUM "Row"

// Flow variables
#define RHO "Density"
#define RHOV "Momentum"
#define RHOE "EnergyStagnationDensity"

#define PS "Pressure"
#define PTREL "Ptrel"
#define PTABS "Ptabs"
#define GENERIC_PT "Pt"
#define TS "Temperature"
#define TTREL "Ttrel"
#define TTABS "Ttabs"
#define GENERIC_TT "Tt"
#define ENTROPY "S_Sref"
#define VELOCITYABS "V"
#define VELOCITYREL "W"
#define GENERIC_VELOCITY "V"
#define MACHABS "Mv"
#define MACHREL "Mw"
#define GENERIC_MACH "M"
#define ALPHA "alpha"
#define BETA "beta"
#define GENERIC_BETA "beta"
#define PHI "phi"
#define MU "mu"
#define SAVED_COORDS "SavedCoordinates"
//
#define HTABS "EnthalpyStagnation"
#define HTREL "RelativeEnthalpyStagnation"

// Geometry variables
#define RADIUS "R"
#define THETA "Theta"

// Reference values (to compute s-sref)
#define PREF 101325.0
#define TREF 288.15


/* Warnings */
#ifdef PRINT_SUPPORT_WARNING
	#define IZAR_WARNING \
		vtkErrorMacro("Avertissement : ce filtre est experimental. A utiliser en connaissance de cause.");
#else
	#define IZAR_WARNING
#endif


#endif
