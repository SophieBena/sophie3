#
# Makefile for Objective-C code base using GNUstep on an Ubuntu cluster
#

export CC=clang
export CXX=clang++
export CXXFLAGS="-std=c++11"
export RUNTIME_VERSION=gnustep-2.0
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
export LD=/usr/bin/ld.gold
export LDFLAGS="-fuse-ld=/usr/bin/ld.gold -L/usr/local/lib"

GNUSTEP_MAKEFILES = /usr/GNUstep/System/Library/Makefiles

# Include the common variables defined by the Makefile Package
include $(GNUSTEP_MAKEFILES)/common.make

OBJ_SWIM_DEW = \
	../thermoengine/thermoengine/aqueous/born.c \
	../thermoengine/thermoengine/aqueous/duanzhang.c \
	../thermoengine/thermoengine/aqueous/holten.c \
	../thermoengine/thermoengine/aqueous/swim.c \
	../thermoengine/thermoengine/aqueous/wagner.c \
	../thermoengine/thermoengine/aqueous/zhangduan.c

OBJ_STEAM = \
	../src/FreeSteam2.1/b23.c \
	../src/FreeSteam2.1/common.c \
	../src/FreeSteam2.1/region1.c \
	../src/FreeSteam2.1/region4.c \
	../src/FreeSteam2.1/steam_Ts.c \
	../src/FreeSteam2.1/steam_ph.c \
	../src/FreeSteam2.1/steam_pv.c \
	../src/FreeSteam2.1/viscosity.c \
	../src/FreeSteam2.1/backwards.c \
	../src/FreeSteam2.1/derivs.c \
	../src/FreeSteam2.1/region2.c \
	../src/FreeSteam2.1/solver2.c \
	../src/FreeSteam2.1/steam_Tx.c \
	../src/FreeSteam2.1/steam_ps.c \
	../src/FreeSteam2.1/surftens.c \
	../src/FreeSteam2.1/zeroin.c \
	../src/FreeSteam2.1/bounds.c \
	../src/FreeSteam2.1/region3.c \
	../src/FreeSteam2.1/steam.c \
	../src/FreeSteam2.1/steam_pT.c \
	../src/FreeSteam2.1/steam_pu.c \
	../src/FreeSteam2.1/thcond.c

OBJC_FILES = \
	../src/AlloyLiquid.m \
	../src/AlloySolid.m \
	../src/BermanAlbite.m \
	../src/BermanProperties+BermanProperties_BermanProperties_ParameterCalibration.m \
	../src/BermanProperties.m \
	../src/BermanStoichiometricPhases.m \
	../src/BiotiteBerman.m \
	../src/ClinoamphiboleBerman.m \
	../src/CpxBerman.m \
	../src/CristobaliteBerman+CristobaliteBerman_ParameterCalibration.m \
	../src/DEWDielectricConstant.m \
	../src/DEWFluid.m \
	../src/DEWH2O.m \
	../src/DEWspecies.m \
	../src/DoubleMatrix.m \
	../src/DoubleTensor.m \
	../src/DoubleVector.m \
	../src/Equilibrate.m \
	../src/EquilibrateState.m \
	../src/EquilibrateUsingMELTSv102.m \
	../src/EquilibrateUsingMELTSv110.m \
	../src/EquilibrateUsingMELTSv120.m \
	../src/EquilibrateUsingStixrude.m \
	../src/EquilibrateUsingpMELTSv561.m \
	../src/EquilibrateUsingMELTSwithDEW.m \
	../src/EquilibrateUsingMELTSandOnlyDEW.m \
	../src/FeldsparBerman+FeldsparBerman_ParameterCalibration.m \
	../src/FeldsparBerman.m \
	../src/FluidDuan.m \
	../src/GarnetBerman.m \
	../src/GehleniteBerman+GehleniteBerman_ParameterCalibration.m \
	../src/GenericH2O.m \
	../src/HKFspeciesComposite.m \
	../src/HKFspeciesProperties.m \
	../src/HKFspeciesProperties+HKFspeciesProperties_ParameterCalibration.m \
	../src/HollandAndPowellProperties.m \
	../src/HollandAndPowellStoichiometricPhases.m \
	../src/HoltenJPCRD2014.m \
	../src/HornblendeBerman.m \
	../src/IntegerVector.m \
	../src/KalsiliteSSBerman.m \
	../src/LeuciteBerman.m \
	../src/LiquidMelts+LiquidMelts_LiquidMelts_ParameterCalibration.m \
	../src/LiquidMelts.m \
	../src/LiquidMeltsCO2+LiquidMeltsCO2_ParameterCalibration.m \
	../src/LiquidMeltsCO2.m \
	../src/LiquidMeltsGenericEM+LiquidMeltsGenericEM_ParameterCalibration.m \
	../src/LiquidMeltsGenericEM.m \
	../src/LiquidMeltsH2O.m \
	../src/LiquidMeltsH2ORevised+LiquidMeltsH2ORevised_ParameterCalibration.m \
	../src/LiquidMeltsH2ORevised.m \
	../src/LiquidMeltsPlusCO2+LiquidMeltsPlusCO2_LiquidMeltsPlusCO2_ParameterCalibration.m \
	../src/LiquidMeltsPlusCO2.m \
	../src/LiquidMeltsPlusOldH2OandNewCO2.m \
	../src/LiquidMeltsSiO2+LiquidMeltsSiO2_ParameterCalibration.m \
	../src/LiquidMeltsSiO2.m \
	../src/LiquidpMelts.m \
	../src/LiquidpMeltsGenericEM.m \
	../src/LiquidpMeltsH2O.m \
	../src/MathSupport.m \
	../src/MeliliteBerman.m \
	../src/Metals.m \
	../src/NephelineSSBerman.m \
	../src/OlivineBerman.m \
	../src/OpxBerman.m \
	../src/OrthoOxideBerman.m \
	../src/OrthoamphiboleBerman.m \
	../src/PhaseBase.m \
	../src/PseudoPhase.m \
	../src/QuartzBerman+QuartzBerman_ParameterCalibration.m \
	../src/RhombohedralBerman.m \
	../src/SanidineBerman+SanidineBerman_ParameterCalibration.m \
	../src/SpinelBerman.m \
	../src/Stixrude.m \
	../src/StixrudeEndmembers.m \
	../src/StixrudeProperties.m \
	../src/StixrudeSolutionPhase.m \
	../src/StixrudeSolutions.m \
	../src/StixrudeStoichiometricPhases.m \
	../src/TridymiteBerman+TridymiteBerman_ParameterCalibration.m \
	../src/Wagner2002.m \
	../src/ZhangDuan2009.m

# Two libraries
LIBRARY_NAME = libswimdew libphaseobjc

#
# The files for the libswimdew library
#
# The C source files to be compiled
libswimdew_C_FILES = $(OBJ_SWIM_DEW) $(OBJ_STEAM)

# The header files for the library
libswimdew_INCLUDE_DIRS = -I/usr/local/include -I../thermoengine/thermoengine/aqueous -I../src/FreeSteam2.1
libswimdew_CFLAGS = -DBUILD_WITH_INLINE_FREESTEAM
libswimdew_LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas

#
# The files for the libphaseobjc library
#
# The C source files to be compiled
libphaseobjc_C_FILES = $(OBJ_STEAM)
# The Objective-C files to compile
libphaseobjc_OBJC_FILES = $(OBJC_FILES)

# The header files for the library
libphaseobjc_INCLUDE_DIRS = -I/usr/local/include -I../src -I../src/FreeSteam2.1
libphaseobjc_CFLAGS = -DBUILD_WITH_INLINE_FREESTEAM
libphaseobjc_OBJCFLAGS = -fobjc-arc -DBUILD_WITH_INLINE_FREESTEAM
libphaseobjc_LDFLAGS = -lgnustep-base -lobjc -L/usr/local/lib -lgsl -lgslcblas -llapack

# Option include to set any additional variables
-include GNUmakefile.preamble

# Include in the rules for making libraries
include $(GNUSTEP_MAKEFILES)/library.make

# Option include to define any additional rules
-include GNUmakefile.postamble
