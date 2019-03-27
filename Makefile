CFLAGS=-std=c89 -pedantic -Wall -Wextra -Wconversion -ftrapv -Wfloat-equal -g -march=native -O3 $(EXTRA_CFLAGS)
LDFLAGS=-g -rdynamic $(EXTRA_LDFLAGS)
LDLIBS=$(EXTRA_LDLIBS)

BIN=denoise

.PHONY: all clean distclean

all: $(BIN)

denoise: denoise.o frame.o dwt.o

denoise.o: denoise.c utils.h config.h

frame.o: frame.c frame.h utils.h

dwt.o: dwt.c frame.h common.h utils.h config.h

clean:
	$(RM) -- *.o $(BIN)

distclean: clean
	$(RM) -- *.gcda
