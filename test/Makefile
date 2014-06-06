
INC_DIR = ../
BUILD_DIR = ../build
OBJS = $(wildcard $(BUILD_DIR)/*)
SRC = $(wildcard test*cpp)
TESTS = $(patsubst %.cpp, u%, $(SRC))
EXEC = $(patsubst %.cpp, %, $(SRC))

CXX = g++
CXXFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lgmp -lm
.PHONY: clean

clean:
	rm -f $(EXEC)

all:	$(TESTS)

utest_%:	$(OBJS)
	if test -e test_$*.cpp; then $(CXX) -o test_$* test_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi
