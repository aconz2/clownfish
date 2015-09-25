CC = g++
CXXFLAGS = -Wall -Werror -std=c++11 -O3 $(shell pkg-config --cflags jellyfish-2.0)
LDFLAGS = $(shell pkg-config --libs-only-L jellyfish-2.0)
LDLIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options $(shell pkg-config --libs-only-l jellyfish-2.0) 

.PHONY: clean all

all: count

count: count.cc 

clean:
	rm -f count

