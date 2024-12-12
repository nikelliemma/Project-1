
# Define variables for object files, source files, headers, output executable, compiler, and flags
SOURCE_PROJECT  = RunTest.cpp Vamana.cpp Graph.cpp dataset.cpp filtered_dataset.cpp
SOURCE_TEST     = filtered_dataset.cpp dataset.cpp Graph.cpp Vamana.cpp unit_testing.cpp 
HEADER          = filtered_dataset.h dataset.h Graph.h Vamana.h acutest.h
OUT_PROJECT     = test_vam
OUT_TEST        = unit_testing
CC              = g++
FLAGS           = -std=c++17 
FLAGS_1         = -std=c++17 -O3 

# Default target to build the executable
all: $(OUT_PROJECT) $(OUT_TEST)

project: $(OUT_PROJECT)

test: $(OUT_TEST)


# Link all source files directly to create the final executable
$(OUT_PROJECT): $(SOURCE_PROJECT) $(HEADER)
	$(CC) $(FLAGS_1) $(SOURCE_PROJECT) -o $(OUT_PROJECT)

#RULE FOR UNIT TESTING
$(OUT_TEST): $(SOURCE_TEST) $(HEADER)
	$(CC) $(FLAGS) $(SOURCE_TEST) -o $(OUT_TEST)


# Clean up generated files
clean:
	rm -f *.o $(OUT_PROJECT) $(OUT_TEST)
