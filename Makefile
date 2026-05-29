.PHONY: all check clean

CXXFLAGS+=-g

all: check

check:
	${MAKE} -C tests check
clean:
	${MAKE} -C tests clean
