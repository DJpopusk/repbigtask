# - inc/
#     - *.h
# - src/
#     - *.c
#     - *.cpp
# - obj/
#     - *.o
# - main

TARGET := main
SOURCES := $(wildcard *.c *.cpp)
OBJECTS := $(patsubst %.c,%.o, $(patsubst %.cpp,%.o,$(SOURCES)))

INCLUDE := -I.
LIBPATH :=
LIBS :=

FLAGS := 
CCFLAGS := $(FLAGS) -std=c99
CXXFLAGS := $(FLAGS)

CC := gcc
Cxx := g++

all: $(OBJECTS)
	$(CC) $(CCFLAGS) $(INCLUDE) $(OBJECTS) -o $(TARGET) $(LIBPATH) $(LIBS)

%.o: %.c
	$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean help

clean:
	rm -rf obj/*
	rm -f $(TARGET)
	
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'