CC = g++
#CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++0x
# remove -Werror for weird warnings from boost 
CXXFLAGS = -Wall  -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++0x
LDFLAGS = $(shell pkg-config --libs jellyfish-2.0) $(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L/-Wl,-rpath=/g') 
LDLIBS=$(shell ls /usr/local/stow/boost-1.54.0/lib/*.so | grep -v python)

all: count
count: count.cc
clean:
	rm -f *.o count
