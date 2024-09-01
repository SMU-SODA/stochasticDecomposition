CPLEXDIR = /Applications/CPLEX_Studio2211/cplex
CONCERTDIR = /Applications/CPLEX_Studio2211/concert
SPALGODIR = /Users/sultan/S6/spAlgorithms

CC = gcc
CFLAGS = -g -O0 -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -ItwoSD_src -I$(SPALGODIR)/solverUtilities -I$(SPALGODIR)/smpsReader -I$(SPALGODIR)/smpsWriter
LDFLAGS = -L$(CPLEXDIR)/lib/arm64_osx/static_pic -L$(CONCERTDIR)/lib/arm64_osx/static_pic -lcplex -lilocplex -lconcert -lm

TARGET = sultan
TEST_TARGET = algo

# Source files, excluding problematic ones
SOURCES = $(wildcard $(SPALGODIR)/solverUtilities/*.c) $(wildcard $(SPALGODIR)/smpsReader/*.c) $(wildcard twoSD_src/*.c)

# Exclude specific files from compilation
EXCLUDE_FILES = $(SPALGODIR)/solverUtilities/solverCheck.c $(SPALGODIR)/smpsReader/readFIles.c twoSD_src/test_algo.c
SOURCES := $(filter-out $(EXCLUDE_FILES), $(SOURCES))

OBJECTS = $(SOURCES:.c=.o)
TEST_SOURCES = twoSD_src/algo.c
TEST_OBJECTS = $(TEST_SOURCES:.c=.o)

all: $(TARGET)

# Linking all object files to create the final executable
$(TARGET): $(OBJECTS)
	$(CC) -w -o $(TARGET) $(OBJECTS) $(LDFLAGS)

# Test target compilation
$(TEST_TARGET): $(OBJECTS) $(TEST_OBJECTS)
	$(CC) -w -o $(TEST_TARGET) $(OBJECTS) $(TEST_OBJECTS) $(LDFLAGS)

# Rule for compiling .c files into .o object files
%.o: %.c
	$(CC) -w $(CFLAGS) -c $< -o $@

# Clean up compiled files
clean:
	rm -f $(TARGET) $(TEST_TARGET) $(OBJECTS) $(TEST_OBJECTS)

# Running the target with the required options
run:
	./sultan -p 20 -i /Users/sultan/S6/spAlgorithms/spInput/ -o /Users/sultan/S6/spAlgorithms/spOutput -e 0 -d 1 -t l -m 1 -c 0
