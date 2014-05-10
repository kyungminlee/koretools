CXXFLAGS = -g -O0 -Wall -std=c++11 -isystem ../include -isystem $(HOME)/.brew/include
LDFLAGS = -L$(HOME)/.brew/lib -lboost_unit_test_framework

bitbox_test: bitbox.cc
	clang++ $(CXXFLAGS) $(LDFLAGS) -o bitbox_test bitbox.cc 
