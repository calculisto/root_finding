DEPENDENCIES_HEADERS=\
	git_submodules/calculisto/array/include \
	git_submodules/calculisto/auto_diff/include \

PROJECT=root_finding
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++2a -Wall -Wextra $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})
LDLIBS+= -lfmt

