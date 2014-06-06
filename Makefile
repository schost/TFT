
INC_DIR = .
BUILD_DIR = build
SOURCES = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

CXX = g++
CXXFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lgmp -lm
.PHONY: clean

all:	$(OBJS)

clean:
	rm -f $(OBJS)
	cd test ; make clean

$(BUILD_DIR)/%.o:	src/%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

