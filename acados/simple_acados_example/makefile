# Simple mass model
# Something wrong with LIBS, check it later
LIBS = 
LIBS += -lacados_c
LIBS += -lblasfeo -lm -lblas -llapack
LIBS += -lacados

MA_OBJS = 
MA_OBJS += mass_model/dyn_mass.o
MA_OBJS += mass_model/power_limit.o

simple_mass_model: $(MA_OBJS) simple_mass_model.cpp
	g++ -o simple_mass_model.out -std=c++11 simple_mass_model.cpp $(MA_OBJS) $(LIBS) -I/usr/local/include/blasfeo/include -L/home/xy-fu/MasterThesis/acados/lib
	@echo  
	@echo "simple mass model build completed"
	@echo  

