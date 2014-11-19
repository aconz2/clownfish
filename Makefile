CC = g++
#CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++0x
# remove -Werror for weird warnings from boost 
CXXFLAGS = -Wall -Werror -Wunused-local-typedefs -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11
# this should get changed for portability sometime!
LDFLAGS = $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib 
LDLIBS = $(shell pkg-config --libs-only-l jellyfish-2.0) -lboost_timer -lboost_chrono -lboost_system -lboost_program_options

all: count
count: count.cc JellyfishHash.hpp
clean:
	rm -f *.o count
