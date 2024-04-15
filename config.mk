DEPENDENCIES_HEADERS= \

PROJECT=root_finding
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++2a -Wall -Wextra $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})
LDLIBS+= -lfmt

