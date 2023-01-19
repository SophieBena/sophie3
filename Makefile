UNAME := $(shell uname)

OBJ_FILES = \
	./src/AlloyLiquid.o \
	./src/AlloySolid.o \
	./src/BermanAlbite.o \
	./src/BermanProperties+BermanProperties_BermanProperties_ParameterCalibration.o \
	./src/BermanProperties.o \
	./src/BermanStoichiometricPhases.o \
	./src/BiotiteBerman.o \
	./src/ClinoamphiboleBerman.o \
	./src/CpxBerman.o \
	./src/CristobaliteBerman+CristobaliteBerman_ParameterCalibration.o \
	./src/DEWDielectricConstant.o \
	./src/DEWFluid.o \
	./src/DEWH2O.o \
	./src/DEWspecies.o \
	./src/DoubleMatrix.o \
	./src/DoubleTensor.o \
	./src/DoubleVector.o \
	./src/Equilibrate.o \
	./src/EquilibrateState.o \
	./src/EquilibrateUsingMELTSv102.o \
	./src/EquilibrateUsingMELTSv110.o \
	./src/EquilibrateUsingMELTSv120.o \
	./src/EquilibrateUsingStixrude.o \
	./src/EquilibrateUsingpMELTSv561.o \
	./src/EquilibrateUsingMELTSwithDEW.o \
	./src/EquilibrateUsingMELTSandOnlyDEW.o \
	./src/FeldsparBerman+FeldsparBerman_ParameterCalibration.o \
	./src/FeldsparBerman.o \
	./src/FluidDuan.o \
	./src/GarnetBerman.o \
	./src/GehleniteBerman+GehleniteBerman_ParameterCalibration.o \
	./src/GenericH2O.o \
	./src/HKFspeciesComposite.o \
	./src/HKFspeciesProperties.o \
	./src/HKFspeciesProperties+HKFspeciesProperties_ParameterCalibration.o \
	./src/HollandAndPowellProperties.o \
	./src/HollandAndPowellStoichiometricPhases.o \
	./src/HoltenJPCRD2014.o \
	./src/HornblendeBerman.o \
	./src/IntegerVector.o \
	./src/KalsiliteSSBerman.o \
	./src/LeuciteBerman.o \
	./src/LiquidMelts+LiquidMelts_LiquidMelts_ParameterCalibration.o \
	./src/LiquidMelts.o \
	./src/LiquidMeltsCO2+LiquidMeltsCO2_ParameterCalibration.o \
	./src/LiquidMeltsCO2.o \
	./src/LiquidMeltsGenericEM+LiquidMeltsGenericEM_ParameterCalibration.o \
	./src/LiquidMeltsGenericEM.o \
	./src/LiquidMeltsH2O.o \
	./src/LiquidMeltsH2ORevised+LiquidMeltsH2ORevised_ParameterCalibration.o \
	./src/LiquidMeltsH2ORevised.o \
	./src/LiquidMeltsPlusCO2+LiquidMeltsPlusCO2_LiquidMeltsPlusCO2_ParameterCalibration.o \
	./src/LiquidMeltsPlusCO2.o \
	./src/LiquidMeltsPlusOldH2OandNewCO2.o \
	./src/LiquidMeltsSiO2+LiquidMeltsSiO2_ParameterCalibration.o \
	./src/LiquidMeltsSiO2.o \
	./src/LiquidpMelts.o \
	./src/LiquidpMeltsGenericEM.o \
	./src/LiquidpMeltsH2O.o \
	./src/MathSupport.o \
	./src/MeliliteBerman.o \
	./src/Metals.o \
	./src/NephelineSSBerman.o \
	./src/OlivineBerman.o \
	./src/OpxBerman.o \
	./src/OrthoOxideBerman.o \
	./src/OrthoamphiboleBerman.o \
	./src/PhaseBase.o \
	./src/PseudoPhase.o \
	./src/QuartzBerman+QuartzBerman_ParameterCalibration.o \
	./src/RhombohedralBerman.o \
	./src/SanidineBerman+SanidineBerman_ParameterCalibration.o \
	./src/SpinelBerman.o \
	./src/Stixrude.o \
	./src/StixrudeEndmembers.o \
	./src/StixrudeProperties.o \
	./src/StixrudeSolutionPhase.o \
	./src/StixrudeSolutions.o \
	./src/StixrudeStoichiometricPhases.o \
	./src/TridymiteBerman+TridymiteBerman_ParameterCalibration.o \
	./src/Wagner2002.o \
	./src/ZhangDuan2009.o

OBJ_SWIM_DEW = \
	./thermoengine/thermoengine/aqueous/born.o \
	./thermoengine/thermoengine/aqueous/duanzhang.o \
	./thermoengine/thermoengine/aqueous/holten.o \
	./thermoengine/thermoengine/aqueous/swim.o \
	./thermoengine/thermoengine/aqueous/wagner.o \
	./thermoengine/thermoengine/aqueous/zhangduan.o

OBJ_SPECIATION = \
	./thermoengine/thermoengine/speciation/nnls.o

OBJ_STEAM = \
	./src/FreeSteam2.1/b23.o \
	./src/FreeSteam2.1/common.o \
	./src/FreeSteam2.1/region1.o \
	./src/FreeSteam2.1/region4.o \
	./src/FreeSteam2.1/steam_Ts.o \
	./src/FreeSteam2.1/steam_ph.o \
	./src/FreeSteam2.1/steam_pv.o \
	./src/FreeSteam2.1/viscosity.o \
	./src/FreeSteam2.1/backwards.o \
	./src/FreeSteam2.1/derivs.o \
	./src/FreeSteam2.1/region2.o \
	./src/FreeSteam2.1/solver2.o \
	./src/FreeSteam2.1/steam_Tx.o \
	./src/FreeSteam2.1/steam_ps.o \
	./src/FreeSteam2.1/surftens.o \
	./src/FreeSteam2.1/zeroin.o \
	./src/FreeSteam2.1/bounds.o \
	./src/FreeSteam2.1/region3.o \
	./src/FreeSteam2.1/steam.o \
	./src/FreeSteam2.1/steam_pT.o \
	./src/FreeSteam2.1/steam_pu.o \
	./src/FreeSteam2.1/thcond.o

ifeq ($(UNAME), Darwin)
all: libphaseobjc.dylib libswimdew.dylib libspeciation.dylib

libspeciation.dylib: $(OBJ_SPECIATION)
	clang -dynamiclib $(OBJ_SPECIATION) -o ./src/libspeciation.dylib

libswimdew.dylib: $(OBJ_SWIM_DEW) $(OBJ_STEAM)
	clang -dynamiclib $(OBJ_SWIM_DEW) $(OBJ_STEAM) -L/usr/local/lib -lgsl -o ./src/libswimdew.dylib

libphaseobjc.dylib: $(OBJ_FILES) $(OBJ_STEAM)
	clang -dynamiclib $(OBJ_FILES) $(OBJ_STEAM) -L/usr/local/lib -lgsl -framework Accelerate -fobjc-arc -fobjc-link-runtime -o ./src/libphaseobjc.dylib

install: libphaseobjc.dylib libswimdew.dylib libspeciation.dylib
	cd src; cp libphaseobjc.dylib /usr/local/lib
	cd src; cp libswimdew.dylib /usr/local/lib
	cd src; cp libspeciation.dylib /usr/local/lib

pyinstall: libswimdew.dylib libspeciation.dylib
	# Need to use special method for conda commands inside makefile
	# conda activate thermoengine
	# conda env update --file environment.yml  --prune
	cd src; cp libswimdew.dylib /usr/local/lib
	cd src; cp libspeciation.dylib /usr/local/lib
	cd thermoengine; pip install --upgrade --use-feature=in-tree-build .

notebook-docs: $(NOTEBOOKS)
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-1-Berman-std-state.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-2-Berman-plus-BM.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-3-HKF.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-4-Stixrude-Debye.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-5-HDNB.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-6-HDNB-Quartz.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-7-Simple-Solution.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-8-Importing-modules.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-9-Complex-Solution.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-10-gas-Speciation-Solution.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-12-SWIM.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-13-Sulfide_liquid.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEW-QuartzSolubility.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEW-Standard-State.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEWFluid.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.0.2-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.0.2-fractionation.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.1.0-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.2.0-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/pMELTS-v5.6.1-adiabatic.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/pMELTS-v5.6.1-melting.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Compare-Phases.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Forsterite-Stixrude.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Phase-Diagram.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Plot-Reaction-Stixrude.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Plot-Reaction.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Quartz-Berman.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Water.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Clinopyroxene-MELTS.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Feldspar-ss-Berman.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Oxide-Geothermometer.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Pyroxene-Geothermometer.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/1-G-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/2-KandC-O2-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/3-K-H2O-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/4-H-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/5-A-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/6-E-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/7a-K-Quartz-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/7b-K-Cor-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/8-K-Qtz-Cor-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/9-K-Qtz-Fld-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/9-Olivine-loop.ipynb

notebook-public: $(NOTEBOOKS)
	cp ./Notebooks/README.md /Users/Shared/PublicNotebooks/
	mkdir -p /Users/Shared/PublicNotebooks/Codegen
	cp -R ./Notebooks/Codegen/* /Users/Shared/PublicNotebooks/Codegen/
	mkdir -p /Users/Shared/PublicNotebooks/Development
	cp -R ./Notebooks/Development/* /Users/Shared/PublicNotebooks/Development/
	mkdir -p /Users/Shared/PublicNotebooks/DEW
	cp -R ./Notebooks/DEW/* /Users/Shared/PublicNotebooks/DEW/
	mkdir -p /Users/Shared/PublicNotebooks/Equilibrate
	cp -R ./Notebooks/Equilibrate/* /Users/Shared/PublicNotebooks/Equilibrate/
	mkdir -p /Users/Shared/PublicNotebooks/MELTS-pMELTS
	cp -R ./Notebooks/MELTS-pMELTS/* /Users/Shared/PublicNotebooks/MELTS-pMELTS/
	mkdir -p /Users/Shared/PublicNotebooks/Presentations
	cp -R ./Notebooks/Presentations/* /Users/Shared/PublicNotebooks/Presentations/
	mkdir -p /Users/Shared/PublicNotebooks/Pure-Phases
	cp -R ./Notebooks/Pure-Phases/* /Users/Shared/PublicNotebooks/Pure-Phases/
	mkdir -p /Users/Shared/PublicNotebooks/Solutions
	cp -R ./Notebooks/Solutions/* /Users/Shared/PublicNotebooks/Solutions/
	mkdir -p /Users/Shared/PublicNotebooks/Calibration
	cp -R ./Notebooks/Calibration/* /Users/Shared/PublicNotebooks/Calibration/

rapid-tests:
	cd thermoengine/thermoengine/test; bash run-rapid-tests.sh
	echo 'Coverage report is available at ./thermoengine/thermoengine/htmlcov/index.html'

tests:
	cd Notebooks; bash test-manual.sh
	echo 'Coverage report is available at ./Notebooks/htmlcov/index.html'

clean:
	cd src; rm -rf *.o *.d ./FreeSteam2.1/*.o ./FreeSteam2.1/*.d libphaseobjc.dylib libswimdew.dylib libspeciation.dylib
	rm -rf ./thermoengine/thermoengine/aqueous/*.o ./thermoengine/thermoengine/aqueous/*.d
	rm -rf ./thermoengine/thermoengine/speciation/*.o ./thermoengine/thermoengine/speciation/*.d

%.o : %.m
	clang -x objective-c -fobjc-arc -I. -I./src/FreeSteam2.1 -DBUILD_WITH_INLINE_FREESTEAM -c $< -o $@

%.o : %.c
	clang -I. -I/usr/local/include -I./thermoengine/thermoengine/aqueous -I./thermoengine/thermoengine/speciation -I./src/FreeSteam2.1 -c $< -o $@
endif


ifeq ($(UNAME), Linux)
all: libphaseobjc.so  libswimdew.so libspeciation.so

libspeciation.so: $(OBJ_SPECIATION)
	clang -shared $(OBJ_SPECIATION) -o ./src/libspeciation.so

libswimdew.so: $(OBJ_SWIM_DEW) $(OBJ_STEAM)
	clang -shared $(OBJ_SWIM_DEW) $(OBJ_STEAM) -L/usr/local/lib -lgsl -lgslcblas -o ./src/libswimdew.so

libphaseobjc.so: $(OBJ_FILES) $(OBJ_STEAM)
	clang -shared $(OBJ_FILES) $(OBJ_STEAM) -lgnustep-base -ldispatch -lgsl -lgslcblas -llapack \
	-L/usr/local/lib -lobjc -fobjc-arc -fobjc-link-runtime -o ./src/libphaseobjc.so

install: libphaseobjc.so libswimdew.so libspeciation.so
	cd src; cp libphaseobjc.so /usr/local/lib
	cd src; cp libswimdew.so /usr/local/lib
	cd src; cp libspeciation.so /usr/local/lib

pyinstall: libswimdew.so libspeciation.so
	cd src; cp libswimdew.so /usr/local/lib
	cd src; cp libspeciation.so /usr/local/lib
	cd thermoengine; pip install --upgrade --use-feature=in-tree-build .

notebook-docs: $(NOTEBOOKS)
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-1-Berman-std-state.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-2-Berman-plus-BM.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-3-HKF.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-4-Stixrude-Debye.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-5-HDNB.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-6-HDNB-Quartz.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-7-Simple-Solution.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-8-Importing-modules.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-9-Complex-Solution.ipynb
#	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-10-gas-Speciation-Solution.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-12-SWIM.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Codegen/Example-13-Sulfide_liquid.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEW-QuartzSolubility.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEW-Standard-State.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/DEW/DEWFluid.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.0.2-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.0.2-fractionation.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.1.0-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/MELTS-v1.2.0-equilibrium.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/pMELTS-v5.6.1-adiabatic.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/MELTS-pMELTS/pMELTS-v5.6.1-melting.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Compare-Phases.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Forsterite-Stixrude.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Phase-Diagram.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Plot-Reaction-Stixrude.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Plot-Reaction.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Quartz-Berman.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Pure-Phases/Water.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Clinopyroxene-MELTS.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Feldspar-ss-Berman.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Oxide-Geothermometer.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Solutions/Pyroxene-Geothermometer.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/1-G-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/2-KandC-O2-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/3-K-H2O-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/4-H-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/5-A-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/6-E-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/7a-K-Quartz-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/7b-K-Cor-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/8-K-Qtz-Cor-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/9-K-Qtz-Fld-Equilibrate.ipynb
	jupyter nbconvert --to rst --execute --output-dir ./Documentation/source ./Notebooks/Equilibrate/9-Olivine-loop.ipynb

notebook-public: $(NOTEBOOKS)
	cp ./Notebooks/README.md /Users/Shared/PublicNotebooks/
	mkdir -p /Users/Shared/PublicNotebooks/Codegen
	cp -R ./Notebooks/Codegen/* /Users/Shared/PublicNotebooks/Codegen/
	mkdir -p /Users/Shared/PublicNotebooks/Development
	cp -R ./Notebooks/Development/* /Users/Shared/PublicNotebooks/Development/
	mkdir -p /Users/Shared/PublicNotebooks/DEW
	cp -R ./Notebooks/DEW/* /Users/Shared/PublicNotebooks/DEW/
	mkdir -p /Users/Shared/PublicNotebooks/Equilibrate
	cp -R ./Notebooks/Equilibrate/* /Users/Shared/PublicNotebooks/Equilibrate/
	mkdir -p /Users/Shared/PublicNotebooks/MELTS-pMELTS
	cp -R ./Notebooks/MELTS-pMELTS/* /Users/Shared/PublicNotebooks/MELTS-pMELTS/
	mkdir -p /Users/Shared/PublicNotebooks/Presentations
	cp -R ./Notebooks/Presentations/* /Users/Shared/PublicNotebooks/Presentations/
	mkdir -p /Users/Shared/PublicNotebooks/Pure-Phases
	cp -R ./Notebooks/Pure-Phases/* /Users/Shared/PublicNotebooks/Pure-Phases/
	mkdir -p /Users/Shared/PublicNotebooks/Solutions
	cp -R ./Notebooks/Solutions/* /Users/Shared/PublicNotebooks/Solutions/
	mkdir -p /Users/Shared/PublicNotebooks/Calibration
	cp -R ./Notebooks/Calibration/* /Users/Shared/PublicNotebooks/Calibration/

tests:
	cd Notebooks; bash test-manual.sh
	echo 'Coverage report is available at ./Notebooks/htmlcov/index.html'

clean:
	cd src; rm -rf *.o *.d ./FreeSteam2.1/*.o ./FreeSteam2.1/*.d libphaseobjc.so libswimdew.so libspeciation.so
	rm -rf ./thermoengine/thermoengine/aqueous/*.o ./thermoengine/thermoengine/aqueous/*.d
	rm -rf ./thermoengine/thermoengine/speciation/*.o ./thermoengine/thermoengine/speciation/*.d

%.o : %.m
	clang -x objective-c -fobjc-runtime=gnustep -fblocks -fobjc-arc -fPIC \
	-I/usr/include/GNUstep -I/usr/include/GNUstep/GNUstepBase -I/usr/local/include/GNUstep \
	-I/usr/local/include -I./thermoengine/thermoengine/aqueous -I./src/FreeSteam2.1 -I. \
	-DBUILD_WITH_INLINE_FREESTEAM -c $< -o $@

%.o : %.c
	clang -I. -I/usr/local/include -I./thermoengine/thermoengine/aqueous -I./thermoengine/thermoengine/speciation -I./src/FreeSteam2.1 -fPIC -c $< -o $@
endif
