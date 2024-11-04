
# Define variables for object files, source files, headers, output executable, compiler, and flags
SOURCE  = test_vam.cpp Vamana.cpp Graph.cpp dataset.cpp
HEADER  = Vamana.h Graph.h dataset.h
OUT     = test_vam
CC      = g++
FLAGS   = -g -O3 -std=c++17

# Default target to build the executable
all: $(OUT)

# Link all source files directly to create the final executable
$(OUT): $(SOURCE) $(HEADER)
	$(CC) $(FLAGS) $(SOURCE) -o $(OUT)

# Clean up generated files
clean:
	rm -f *.o $(OUT)
