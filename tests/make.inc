CXX = g++-4.9
INCLUDE_PATHS = ../include $(HOME)/.brew/include
LIBRARY_PATHS = $(HOME)/.brew/lib
LIBRARIES = # boost_unit_test_framework-mt

CXXFLAGS = -g -O0 -Wall -std=c++11
LDFLAGS = $(HOME)/.brew/lib/libboost_unit_test_framework-mt.a
