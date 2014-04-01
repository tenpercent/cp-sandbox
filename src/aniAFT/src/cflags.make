#Makefile inlude
CC  = gcc
CXX = g++
AR  = ar

LDFLAGS  := $(LDFLAGS) -L$(AFTHOME)/lib
LDLIBS   := $(LDLIBS) -lm
FLAGS    := -W -Wall -O3 -march=native
CFLAGS   := $(CFLAGS) $(FLAGS)
CXXFLAGS := $(CXXFLAGS) $(FLAGS)

############################################################
version = $(shell cat $(AFTHOME)/../../VERSION)
AFTBIN  = $(AFTHOME)/../../bin
AFTLIB  = $(AFTHOME)/../../lib
AFTDAT  = $(AFTHOME)/../../data
LIBAFT  = $(AFTLIB)/libaft3D-$(version).a

