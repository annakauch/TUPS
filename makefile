DIR = ./

#compiler
CC = mpif90 

#Flags - ipo takes long
#production flags
FFLAGS = -Ofast -ipo -g -heap-arrays -xHost -traceback -check bounds -module MOD

#developement flags
#FFLAGS = -O1 -heap-arrays -xHost -traceback -check bounds -module MOD

#for debugging
#FFLAGS = -O0 -ggdb -xHost -heap-arrays -traceback -check bounds -warn all -module MOD

#load libraries
LIB = $(DIR)/lib/libdfftpack.a -mkl -lhdf5_fortran -lhdf5hl_fortran -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

OBJ_DIR = OBJ
SRC_DIR = SRC

OBJ_FILES = $(addprefix $(OBJ_DIR)/, hdf5_wrapper.o global_parameter.o MPI_mod.o parquet_ini.o MATH_mod.o parquet_formfactors.o parquet_util.o parquet_plot.o parquet_chi.o parquet_BSE.o parquet_PhiR.o parquet_equation.o parquet_EQMatrix.o parquet_EQContribs.o parquet_EQOuter.o parquet_selfenergy.o parquet_SDE.o parquet_sus.o parquet_sus_opt.o parquet_check.o preprocessing.o loop.o tups.o)
SRC_FILES = $(addprefix $(SRC_DIR)/, hdf5_wrapper.f90 global_parameter.f90 MPI_mod.f90 parquet_ini.f90 MATH_mod.f90 parquet_formfactors.f90 parquet_util.f90 parquet_plot.f90 parquet_chi.f90 parquet_BSE.f90 parquet_PhiR.f90 parquet_equation.f90 parquet_EQMatrix.f90 parquet_EQContribs.f90 parquet_EQOuter.f90 parquet_selfenergy.f90 parquet_SDE.f90 parquet_sus.f90 parquet_sus_opt.f90 parquet_check.f90 preprocessing.f90 loop.o tups.f90)

TUPS: $(OBJ_FILES)
	$(CC) $(FFLAGS) $(OBJ_FILES) $(LIB) -o TUPS

$(OBJ_DIR)/tups.o: $(SRC_DIR)/tups.f90 $(SRC_DIR)/hdf5_wrapper.f90 $(SRC_DIR)/MPI_mod.f90 $(SRC_DIR)/global_parameter.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_plot.f90 $(SRC_DIR)/parquet_chi.f90 $(SRC_DIR)/parquet_BSE.f90 $(SRC_DIR)/parquet_equation.f90 $(SRC_DIR)/parquet_selfenergy.f90 $(SRC_DIR)/parquet_sus.f90 $(SRC_DIR)/parquet_sus_opt.f90  $(SRC_DIR)/parquet_check.f90 $(SRC_DIR)/preprocessing.f90 
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/loop.o: $(SRC_DIR)/loop.f90 $(SRC_DIR)/hdf5_wrapper.f90 $(SRC_DIR)/MPI_mod.f90 $(SRC_DIR)/global_parameter.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_plot.f90 $(SRC_DIR)/parquet_chi.f90 $(SRC_DIR)/parquet_BSE.f90 $(SRC_DIR)/parquet_equation.f90 $(SRC_DIR)/parquet_selfenergy.f90 $(SRC_DIR)/parquet_check.f90 $(SRC_DIR)/preprocessing.f90 
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/preprocessing.o: $(SRC_DIR)/preprocessing.f90 $(SRC_DIR)/parquet_check.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_chi.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_plot.f90 $(SRC_DIR)/hdf5_wrapper.f90 $(SRC_DIR)/global_parameter.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 

$(OBJ_DIR)/parquet_check.o: $(SRC_DIR)/parquet_check.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_ini.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 

$(OBJ_DIR)/parquet_sus.o: $(SRC_DIR)/parquet_sus.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_BSE.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 

$(OBJ_DIR)/parquet_SDE.o: $(SRC_DIR)/parquet_SDE.f90 $(SRC_DIR)/global_parameter.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_BSE.f90 $(SRC_DIR)/parquet_equation.f90
	$(CC) $(FFLAGS)  -c $< -o $@  
 
$(OBJ_DIR)/parquet_selfenergy.o: $(SRC_DIR)/parquet_selfenergy.f90 $(SRC_DIR)/global_parameter.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_BSE.f90 $(SRC_DIR)/parquet_equation.f90
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/parquet_EQOuter.o: $(SRC_DIR)/parquet_EQOuter.f90 $(SRC_DIR)/parquet_EQContribs.f90 $(SRC_DIR)/parquet_EQMatrix.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_BSE.f90
	$(CC) $(FFLAGS) -c $< -o $@  

$(OBJ_DIR)/parquet_EQContribs.o: $(SRC_DIR)/parquet_EQContribs.f90 $(SRC_DIR)/parquet_EQMatrix.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_BSE.f90
	$(CC) $(FFLAGS) -c $< -o $@  

$(OBJ_DIR)/parquet_EQMatrix.o: $(SRC_DIR)/parquet_EQMatrix.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_BSE.f90
	$(CC) $(FFLAGS) -c $< -o $@  

$(OBJ_DIR)/parquet_equation.o: $(SRC_DIR)/parquet_equation.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_BSE.f90
	$(CC) $(FFLAGS) -c $< -o $@  

$(OBJ_DIR)/parquet_PhiR.o: $(SRC_DIR)/parquet_PhiR.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90
	$(CC) $(FFLAGS) -c $< -o $@  

$(OBJ_DIR)/parquet_BSE.o: $(SRC_DIR)/parquet_BSE.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/parquet_util.f90
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/parquet_chi.o: $(SRC_DIR)/parquet_chi.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/global_parameter.f90 $(SRC_DIR)/MATH_mod.f90
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/parquet_plot.o: $(SRC_DIR)/parquet_plot.f90 $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/hdf5_wrapper.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 
 
$(OBJ_DIR)/parquet_util.o: $(SRC_DIR)/parquet_util.f90 $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/global_parameter.f90
	$(CC) $(FFLAGS)  -c $< -o $@  

$(OBJ_DIR)/parquet_formfactors.o: $(SRC_DIR)/parquet_formfactors.f90 $(SRC_DIR)/parquet_ini.f90
	$(CC) $(FFLAGS)  -c $< -o $@

$(OBJ_DIR)/MATH_mod.o: $(SRC_DIR)/MATH_mod.f90 $(SRC_DIR)/parquet_ini.f90 $(SRC_DIR)/global_parameter.f90
	$(CC) $(FFLAGS)  -c $< -o $@
	
$(OBJ_DIR)/parquet_ini.o: $(SRC_DIR)/parquet_ini.f90
	$(CC) $(FFLAGS)  -c $< -o $@
	
$(OBJ_DIR)/MPI_mod.o: $(SRC_DIR)/MPI_mod.f90 $(SRC_DIR)/global_parameter.f90
	$(CC) $(FFLAGS)  -c $< -o $@

$(OBJ_DIR)/global_parameter.o: $(SRC_DIR)/global_parameter.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 

$(OBJ_DIR)/hdf5_wrapper.o: $(SRC_DIR)/hdf5_wrapper.f90
	$(CC) $(FFLAGS)  -c $< -o $@ 

clean: 
	-rm $(OBJ_FILES) MOD/*.mod TUPS
