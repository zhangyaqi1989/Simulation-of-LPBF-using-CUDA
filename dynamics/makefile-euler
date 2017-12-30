# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT		:= -O3
ARCH   	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11

# Linker options
LDOPT 	:= $(OPT)
LDFLAGS := 
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : cpu_dynamic gpu_dynamic


gpu_dynamic: kernel.cu 
	module load cuda;nvcc -o gpu_dynamic $(OPT) kernel.cu -ccbin $(BIN)

cpu_dynamic: deposition_sim.cpp
	g++ -O3 -g -o cpu_dynamic deposition_sim.cpp
# TODO: add targets for building executables

.PHONY: clean
clean:
	rm -f gpu_dynamic cpu_dynamic
