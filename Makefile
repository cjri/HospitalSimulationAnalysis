CC	      = gcc
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
BAS		= process_simulation.o io.o utilities.o

simulation: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_sim_process  $(LD_FLAGS)
process_simulation.o: process_simulation.cpp
	$(CC) $(CC_FLAGS) -c process_simulation.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

