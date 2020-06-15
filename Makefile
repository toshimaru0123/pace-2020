CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11 

PROGEAM = wankoSOBA

all: $(PROGEAM)

$(PROGEAM): exactTDD.cpp
	$(CXX) $(CXXFLAGS) -o $(PROGEAM) exactTDD.cpp