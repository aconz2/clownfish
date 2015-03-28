CC = g++
CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11
# this is very unportable
LDFLAGS := $(LDFLAGS) $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib 
LDLIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options  $(shell pkg-config --libs-only-l jellyfish-2.0) -larmadillo 

all: bin/clownfish-count bin/clownfish-index bin/clownfish-solve
bin/clownfish-count: count
	cp $< $@
bin/clownfish-index: index
	cp $<  $@
bin/clownfish-solve: solve
	cp $< $@

index: index.cc
solve: solve.cc
count: count.cc 
clean:
	rm -f bin/* index solve count
