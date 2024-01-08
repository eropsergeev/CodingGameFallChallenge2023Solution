CXX ?= g++
CXXFLAGS ?= -std=c++23 -DLOCAL -Wall -Wextra -g -Wno-parentheses -pipe -O0 -fsanitize=address -fsanitize=undefined  -mavx2 -mbmi2 -mpopcnt
OPTFLAGS ?= -std=c++23 -pipe -Ofast -fpic -fno-rtti -fno-exceptions -fcf-protection=none -fno-stack-protector -march=native
SUBMISSIONFLAGS ?= -std=c++23 -pipe -Ofast -fpic -fno-rtti -fno-exceptions -fcf-protection=none -fno-stack-protector -mavx2 -mbmi2 -mpopcnt
#  -fallow-store-data-races -ffast-math -funroll-loops -falign-functions -falign-jumps -falign-labels -falign-loops -freorder-blocks-algorithm=stc -fno-unroll-loops

all: main optimized

main: main.cpp $(wildcard *.hpp)
	$(CXX) -o $@ $< $(CXXFLAGS)

submission.S: main.cpp $(wildcard *.hpp)
	$(CXX) -o $@ $< $(SUBMISSIONFLAGS) -S

compress: compress.cpp
	$(CXX) -o $@ $< $(OPTFLAGS)

submission.cpp:	submission.S compress
	./compress < submission.S > submission.cpp

optimized: submission.cpp
	$(CXX) -pipe -std=c++17 -o optimized $<

