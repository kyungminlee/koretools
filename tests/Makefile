include ./make.inc

TARGETS = bitbox_test collector_test

INCLUDE_FLAGS = $(foreach d, $(INCLUDE_PATHS), -isystem$d)
LIBRARY_FLAGS = $(foreach d, $(LIBRARY_PATHS), -L$d) $(foreach l, $(LIBRARIES), -l$l)

all: $(TARGETS)

% : %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LDFLAGS) $(LIBRARY_FLAGS) $< -o $@  
clean :
	rm -rf $(TARGETS) *.dSYM
