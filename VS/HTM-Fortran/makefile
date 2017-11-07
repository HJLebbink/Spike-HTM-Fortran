#filename Spike\Spike-Fortran-Main\makefile

IFORT 			= ifort
GFORTRAN 		= gfortran
AR 				= ar
AR_MIC          = /usr/linux-k1om-4.7/bin/x86_64-k1om-linux-ar
OBJDUMP			= /usr/linux-k1om-4.7/bin/x86_64-k1om-linux-objdump

LOCAL_PATH 		= /home/henk/sc/ScHJ/CWS

INCPATH 		= -I$(LOCAL_PATH)
OBJPATH 		= 
#REPORT 		= -qopt-report=5 -qopt-report-phase=vec -qopt-report-phase=par -qopenmp-report=2 
REPORT			=

MAIN            = main_spike_fortran
SOURCES 		= 	./src/rdtsc.s 					\
					./src/tools_m.f90 				\
					./src/constants_m.f90 			\
					./src/spike_tools_m.f90			\
					./src/time_tracer_m.f90			\
					./src/performance_c.f90			\
					./src/dumper_spike_c.f90 		\
					./src/dumper_state_c.f90		\
					./src/dumper_topology_c.f90		\
					./src/network_c.f90				\
					./src/topology_analysis_m.f90
				  
CFLAGS_INTEL    	 = -fpp -standard-semantics -qopenmp -names as_is -warn all -module ./mod/ $(INCPATH) $(REPORT)
CFLAGS_INTEL_DEBUG	 = -O0 -g -check all -traceback -fpe-all=0 $(CFLAGS_INTEL)

CFLAGS_INTEL_RELEASE = -O1 -xCORE-AVX2 -no-prec-div -fp-model fast=2 $(CFLAGS_INTEL)
#CFLAGS_INTEL_RELEASE = -O3 -xCORE-AVX2 -no-prec-div -fp-model fast=2 $(CFLAGS_INTEL)
#CFLAGS_INTEL_RELEASE = -O3 -ipo -no-prec-div -static -fp-model fast=2 -parallel -guide $(CFLAGS_INTEL)
CFLAGS_MIC_RELEASE   = -mmic -O1 -ipo -no-prec-div -static -fp-model fast=2 $(CFLAGS_INTEL)

DUMP_ASM 			= -ipo-S -use-msasm -fsource-asm -fcode-asm -Fa./asm/ 


CFLAGS_GNU           = -cpp -Wall $(INCPATH)
CFLAGS_GNU_RELEASE   = -O2 $(CFLAGS_GNU)
CFLAGS_GNU_DEBUG     = -Og -g -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace $(CFLAGS_GNU)

#LIBS				 = $(LOCAL_PATH)/Spike/Spike-Fortran-LIB/libspikefortran.a
LIBS			     =

all: intel

asm: generate_asm_intel
mic: clean compile_mic_native copy run
mic2: clean compile_mic_offload
gnu: clean compile_gnu
intel: clean compile_intel
intel_debug: clean compile_intel_debug


generate_asm_intel:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL X64: Generating ASM (RELEASE) --------------------------------"
	@$(IFORT) $(CFLAGS_INTEL_RELEASE) $(DUMP_ASM) $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out
	@$(OBJDUMP) --disassemble --all-headers --demangle=gnu --file-start-context  -Mintel,intel-mnemonic --source $(MAIN).out > ./asm/$(MAIN).asm

compile_intel:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL X64: Compiling + Linking (RELEASE) ---------------------------"
	$(IFORT) $(CFLAGS_INTEL_RELEASE) -qoffload=none -qno-openmp-offload $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out

compile_intel_debug:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL X64: Compiling + Linking (DEBUG) -----------------------------"
	$(IFORT) $(CFLAGS_INTEL_DEBUG) -qoffload=none -qno-openmp-offload $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out

compile_mic_native:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL MIC Native: Compiling + Linking (RELEASE) --------------------"
	$(IFORT) $(CFLAGS_MIC_RELEASE) -qoffload=none -qno-openmp-offload $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out

compile_mic_offload:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL MIC Offload: Compiling + Linking (RELEASE) -------------------"
	$(IFORT) $(CFLAGS_INTEL_RELEASE) -qoffload=mandatory $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out

compile_mic_offload_debug:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL MIC Offload: Compiling + Linking (RELEASE) -------------------"
	$(IFORT) $(CFLAGS_MIC_DEBUG) -qoffload=mandatory $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out

generate_asm_mic:
	@echo "--------------------------------------------------------------------"
	@echo "INTEL X64: Generating MIC (RELEASE) --------------------------------"
	@$(IFORT) $(CFLAGS_INTEL_RELEASE) $(DUMP_ASM) $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out
	@$(OBJDUMP) --disassemble --all-headers --demangle=gnu --file-start-context  -Mintel,intel-mnemonic --source --architecture=k1om $(MAIN).out > ./asm/$(MAIN).asm

	
compile_gnu:
	@echo "--------------------------------------------------------------------"
	@echo "GNU X64: Compiling + Linking (RELEASE) -----------------------------"
	$(GFORTRAN) $(CFLAGS_GNU_RELEASE) $(SOURCES) ./src/$(MAIN).f90 $(LIBS) -o $(MAIN).out
		
copy: 
	@echo "password TamK76gZ"
	@sudo scp  $(MAIN).out mic0:
	@sudo dos2unix run.sh
	@sudo scp run.sh mic0:
	@sudo scp /opt/intel/composerxe/lib/mic/libiomp5.so mic0:/tmp

run: 
	@echo "do manually:"
	@echo "./run.sh"
	@echo "exit"
	sudo ssh mic0
	
clean: 
	@-rm -f ./mod/*.mod
	@-rm -f $(MAIN).out
	@-rm -f *.optrpt
	
