# Warnings
WFLAGS := -Wall

# Optimization
OPT := -O3

# Language standard
CXXSTD := -std=c++11
BIN = "/usr/local/gcc/6.4.0/bin/gcc"

EXEC := main cpu_thermal gpu_thermal

all: $(EXEC) makefile

%.o : %.cpp %.hpp makefile constants.hpp
	@echo Compiling $<...
	$(CXX) $(CXXSTD) $(WFLAGS) $(OPT) -c $< -o $@

main: main.o particle.o dynamics.o thermal.o utility.o constants.hpp
	@echo Compiling executable $@...
	$(CXX) $(CXXSTD) $(WFLAGS) $(OPT) -o $@ $^

cpu_thermal: cpu_thermal.o particle.o dynamics.o thermal.o utility.o
	@echo Compiling executable $@...
	$(CXX) $(CXXSTD) $(WFLAGS) $(OPT) -o $@ $^

cuda_thermal.o: cuda_thermal.cu cuda_thermal.cuh
	@echo Compiling cuda_thermal.o...
	module load cuda; nvcc $(OPT) -c cuda_thermal.cu -ccbin $(BIN)

gpu_thermal.o: gpu_thermal.cu constants.hpp
	@echo Compiling gpu_thermal.o...
	module load cuda; nvcc $(OPT) -c gpu_thermal.cu -ccbin $(BIN)

gpu_thermal: gpu_thermal.o particle.o utility.o dynamics.o
	@echo Compiling executable $@...
	module load cuda; nvcc $(OPT) -o $@ $^ -ccbin $(BIN)

clean:
	rm $(EXEC) *.o
