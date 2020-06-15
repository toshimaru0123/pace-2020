CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11 

all: exact

exact: exactTDD.cpp
	$(CXX) $(CXXFLAGS) -o exactTDD exactTDD.cpp