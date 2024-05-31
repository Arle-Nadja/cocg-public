COCG_k = objects/cocg_kernel.o
COCG_c = objects/cocg_components.o 
COCG_p = objects/cocg_procedures.o
modin = -Imodules 
modout = -Jmodules 
lapack = -LC:\lapack-3.9.1 -llapack -lrefblas 

all : mkdr a.exe

mkdr : 
ifeq ($(OS),Windows_NT)
	@if not exist objects mkdir objects
else
	@mkdir -p objects
endif

ifeq ($(OS),Windows_NT)
	@if not exist modules mkdir modules
else
	@mkdir -p modules
endif

a.exe : CG_test.f90 $(COCG_k) $(COCG_c) $(COCG_p)
	gfortran CG_test.f90 $(COCG_k) $(COCG_c) $(COCG_p) $(modin)

$(COCG_k) : cocg_kernel.f90 $(COCG_p)
	gfortran -c -o $(COCG_k) cocg_kernel.f90 $(modout)

$(COCG_p) : cocg_procedures.f90 $(COCG_c)
	gfortran -c -o $(COCG_p) cocg_procedures.f90 $(modout)

$(COCG_c) : cocg_components.f90
	gfortran -c -o $(COCG_c) cocg_components.f90 $(modout)

clean : 
	del /S *.mod *.o *.exe