CC := gcc
CXX := g++
CFLAGS := -std=c++11 -shared -fPIC
GRAPHLAB_SDK_HOME := ~/workspace/graphlab-sdk

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CFLAGS += --stdlib=libc++ -undefined dynamic_lookup
endif

SRCS := $(wildcard *.cpp)
TARGETS := $(SRCS:%.cpp=%.so)

all: $(TARGETS)

%.so: %.cpp
	$(CXX) -o $@ $(CFLAGS) -I $(GRAPHLAB_SDK_HOME) $^

clean:
	rm *.so
