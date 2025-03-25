###############################################################################
# File: Makefile
# Author: Fatima-Zahra Mihami
# Date: 2023-10-27
###############################################################################

DEBUG ?= 0

INCLUDE_DIR = include
SOURCE_DIR = src
BUILD_DIR = build

#####################
## COMPILER CALLS  ##
#####################

CXX = g++
CXXFLAGS = -std=c++11 -fopenmp
CPPFLAGS = -I $(INCLUDE_DIR)
LDFLAGS = -fopenmp

# Debug Flags
ifeq ($(DEBUG), 1)
CXXFLAGS += -g -O0 -Wall -Wextra -pedantic
EXECUTABLE_NAME = barracuda_debug
else
CXXFLAGS += -O3 -w
EXECUTABLE_NAME = barracuda
endif

CXX_COMPILER_CALL = $(CXX) $(CXXFLAGS) $(CPPFLAGS)

CXX_SOURCES = $(wildcard $(SOURCE_DIR)/*.cpp)
CXX_OBJECTS = $(patsubst $(SOURCE_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(CXX_SOURCES))

##############
## TARGETS  ##
##############
all: create build

create:
	@mkdir -p build

build: $(CXX_OBJECTS)
	$(CXX_COMPILER_CALL) $(CXX_OBJECTS) $(LDFLAGS) -o $(EXECUTABLE_NAME)

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(EXECUTABLE_NAME)


##################
## COMPILE CPP  ##
##################
$(BUILD_DIR)/main.o: $(SOURCE_DIR)/main.cpp
	$(CXX_COMPILER_CALL) -c $(SOURCE_DIR)/main.cpp -o $(BUILD_DIR)/main.o

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(INCLUDE_DIR)/%.h
	$(CXX_COMPILER_CALL) -c $< -o $@


###########
## PHONY ##
###########
.PHONY: create build clean