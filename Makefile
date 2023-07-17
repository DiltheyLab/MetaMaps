CXXFLAGS += -O2 -ggdb -fopenmp -std=c++11 -Isrc -L /usr/lib/x86_64-linux-gnu//lib -I /usr/lib/x86_64-linux-gnu//include 
CPPFLAGS += -DUSE_BOOST

UNAME_S=$(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CXXFLAGS += -mmacosx-version-min=10.7 -stdlib=libc++
else
	CXXFLAGS += -include src/common/memcpyLink.h -Wl,--wrap=memcpy
	CFLAGS += -include src/common/memcpyLink.h
endif

SOURCES=src/map/mash_map.cpp

OBJECTS=$(SOURCES:.cpp=.o)  

all : metamaps

metamaps : $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OBJECTS) -o metamaps -lstdc++ -fopenmp -lz -lm -lpthread -lboost_system -lboost_filesystem -lboost_serialization -lboost_regex -lboost_math_c99

.SUFFIXES :

%.o : %.cpp 
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o : %.c++
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

install : metamaps
	mkdir -p /usr/local/bin/
	cp `pwd`/metamaps /usr/local/bin/

clean :
	-rm -f metamaps
	-rm -f src/map/*.o
