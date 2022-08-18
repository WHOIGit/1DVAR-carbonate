# Makefile for WAP-1D-VAR v1.0

MAIN=	driver.o
OBJS= 	common_mod.o		cost.o		\
	derivs_mod.o 		eco_common.o		eco_derivs.o	\
	eco_params.o		forcing.o		grid.o		\
	io.o					light.o		\
	model_mod.o		numeric_subs.o		physderivs_mod.o\
	types.o
ADJMAIN=adjoint_driver.o
ADJOBS=	adcost.o		adderivs_mod.o		adeco_common.o	\
	adeco_derivs.o		adj_common.o				\
	adlight.o		admodel_mod.o		adnumeric_subs.o\
	adphysderivs_mod.o	hessian_mod.o		m1qn3_mod.o	\
	adBuffer.o		adStack.o
ADJTEST=adjoint_test.o
HESSIAN=hessian_driver.o
OPTMAIN=	opt_driver.o
OPTOBJS= 	common_mod.o		cost.o		\
	derivs_mod.o 		eco_common.o		eco_derivs.o	\
	eco_params.o		forcing.o		grid.o		\
	io.o					light.o		\
	model_mod.o		numeric_subs.o		physderivs_mod.o\
	types.o		newuoa_mod.o		opt_common.o

# Fortran compiler
FC = mpif90
# standard optimizations, no debugging:
FFLAGS = -O3 -fast -r8 $(DFLAGS)
FFLAGS1 = -r8 $(DFLAGS)
# need to point to netcdf.mod
INCLUDES = -I/vortexfs1/apps/pgi/linux86-64/2018/netcdf/netcdf-4.5.0/openmpi-2.1.2/include
LIBNETCDF = /vortexfs1/apps/pgi/linux86-64/2018/netcdf/netcdf-4.5.0/openmpi-2.1.2/lib
LIBHDF = /vortexfs1/apps/pgi/linux86-64/2018/netcdf/netcdf-4.5.0/openmpi-2.1.2/lib
LIBS = -L$(LIBNETCDF) -L$(LIBHDF) -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5 
RM = rm -f

VPATH = ../framework

# compile rules
.SUFFIXES:
.SUFFIXES: .F .F90 .f90 .f .for .c .o


.F.o: 
	$(FC) $(FFLAGS) $(INCLUDES) -c   $<
.F90.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c   $<
.f.o: 
	$(FC) $(FFLAGS) $(INCLUDES) -c   $<
.for.o: 
	$(FC) $(FFLAGS) $(INCLUDES) -c   $<
.f90.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c   $<
.c.o:
	$(CC) -c   $<
.f:
	$(FC)  $(FFLAGS) $(INCLUDES) $< -o  $@
.for:
	$(FC)  $(FFLAGS) $(INCLUDES) $< -o  $@
.f90:
	$(FC)  $(FFLAGS) $(INCLUDES) $< -o  $@
.c:
	$(CC)  $< -o  $@

all: driver adjoint_driver test_adjoint hessian_driver opt_driver

driver: $(MAIN) $(OBJS)
	$(FC) $(FFLAGS1) -o driver $(MAIN) $(OBJS) $(LIBS)

adjoint_driver: $(ADJMAIN) $(OBJS) $(ADJOBS)
	$(FC) $(FFLAGS1) -o adjoint_driver $(ADJMAIN) $(OBJS) $(ADJOBS) $(LIBS)

test_adjoint: $(ADJTEST) $(OBJS) $(ADJOBS)
	$(FC) $(FFLAGS1) -o test_adjoint $(ADJTEST) $(OBJS) $(ADJOBS) $(LIBS)

hessian_driver: $(HESSIAN) $(OBJS) $(ADJOBS)
	$(FC) $(FFLAGS1) -o hessian_driver $(HESSIAN) $(OBJS) $(ADJOBS) $(LIBS)

opt_driver: $(OPTMAIN) $(OPTOBJS)
	$(FC) $(FFLAGS1) -o opt_driver $(OPTMAIN) $(OPTOBJS) $(LIBS)

adcost.o:		const.o			cost.o		\
			grid.o
adderivs_mod.o:		adeco_common.o		const.o		\
			eco_common.o		eco_params.o	\
			adBuffer.o		adStack.o
adeco_common.o:		grid.o			
adeco_derivs.o:		adderivs_mod.o		adeco_common.o	\
			adlight.o				\
			const.o			cost.o		\
			eco_common.o		eco_params.o	\
			forcing.o		grid.o
adj_common.o:
adjoint_driver.o: 	adj_common.o		admodel_mod.o	\
			common_mod.o		const.o		\
			cost.o			eco_common.o	\
			eco_derivs.o		eco_params.o	\
			forcing.o		grid.o		\
			hessian_mod.o		io.o		\
			light.o		\
			model_mod.o		m1qn3_mod.o		
adjoint_test.o: 	adj_common.o		admodel_mod.o	\
			common_mod.o		const.o		\
			cost.o			eco_common.o	\
			eco_derivs.o		eco_params.o	\
			forcing.o		grid.o		\
			io.o			light.o		\
			model_mod.o		m1qn3_mod.o
adlight.o:		const.o			grid.o		\
			light.o
admodel_mod.o:		adcost.o		adeco_derivs.o	\
			adphysderivs_mod.o	common_mod.o	\
			const.o			cost.o		\
			eco_derivs.o		eco_params.o	\
			forcing.o		grid.o		\
			io.o			model_mod.o	\
			physderivs_mod.o
adnumeric_subs.o:	const.o	
adphysderivs_mod.o:	adnumeric_subs.o	common_mod.o	\
			const.o			eco_common.o	\
			eco_params.o		grid.o		\
			physderivs_mod.o
common_mod.o:		const.o			eco_params.o	
const.o:		
cost.o:			common_mod.o		const.o		\
			forcing.o		grid.o	
derivs_mod.o:		const.o			eco_common.o	\
			eco_params.o			io.o
driver.o: 		common_mod.o		cost.o		\
			eco_common.o		eco_derivs.o	\
			eco_params.o		forcing.o	\
			grid.o			io.o		\
			light.o			model_mod.o
eco_common.o:		common_mod.o		cost.o		\
			eco_params.o		grid.o		\
					types.o
eco_derivs.o:	 	common_mod.o		const.o		\
			cost.o			derivs_mod.o	\
			eco_common.o		eco_params.o	\
			forcing.o		grid.o		\
			light.o
eco_params.o:		
forcing.o:		common_mod.o		const.o		\
			grid.o
grid.o:			common_mod.o		const.o		\

hessian_driver.o:	admodel_mod.o		common_mod.o	\
			const.o		\
			cost.o			eco_common.o	\
			eco_derivs.o		eco_params.o	\
			forcing.o		grid.o		\
			io.o			light.o		\
			model_mod.o
hessian_mod.o:		const.o			eco_params.o	\

io.o:			common_mod.o		cost.o		\
			eco_common.o		eco_params.o	\
			forcing.o		grid.o		\
			types.o
light.o: 		const.o			grid.o		\

model_mod.o: 		common_mod.o		const.o		\
			cost.o			eco_derivs.o	\
			eco_params.o		forcing.o	\
			grid.o			io.o		\
			physderivs_mod.o
m1qn3_mod.o:		hessian_mod.o
numeric_subs.o:		
physderivs_mod.o:	common_mod.o		const.o		\
			eco_common.o		eco_params.o	\
			forcing.o		grid.o		\
			numeric_subs.o
types.o:		
adBuffer.o:
adStack.o:
opt_driver.o: 		common_mod.o		cost.o		\
			eco_common.o		eco_derivs.o	\
			eco_params.o		forcing.o	\
			grid.o			io.o		\
			light.o			model_mod.o	\
			newuoa_mod.o		opt_common.o
newuoa_mod.o:	model_mod.o	eco_params.o	io.o
opt_common.o: 		common_mod.o	eco_params.o
clean: 
	$(RM) *.o *.oo *.mod driver adjoint_driver test_adjoint hessian_driver opt_driver
