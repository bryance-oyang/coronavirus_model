EXEC=z.out
srcdir=

SHELL=/bin/sh
CC=gcc -pipe -mtune=native -march=native
OFLAGS=-O3 -flto -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE
CFLAGS+=-std=gnu90 -Wall -Wextra
LDFLAGS=-lc -lm
CDEBUG=-g -p -DCDEBUG
DFLAGS=$(CFLAGS) -M

ifdef srcdir
VPATH=$(srcdir)
SRCS=$(wildcard $(srcdir)/*.c)
HDRS=$(wildcard $(srcdir)/*.h)
CFLAGS+=-I. -I$(srcdir)
else
SRCS=$(wildcard *.c)
HDRS=$(wildcard *.h)
endif
OBJS=$(SRCS:.c=.o)
DEPS=$(SRCS:.c=.d)
ASMS=$(SRCS:.c=.s)

ifeq ($(MAKECMDGOALS), debug)
CFLAGS+=$(CDEBUG)
else
CFLAGS+=$(OFLAGS)
LDFLAGS+=$(OFLAGS)
endif

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

.DEFAULT_GOAL=all
.PHONY: all
all: $(DEPS) $(EXEC)
	@echo done

.PRECIOUS: data/infected_00000.dat
data/infected_00000.dat: $(DEPS) $(EXEC)
	-rm -rf data
	mkdir data
	-./z.out

.PRECIOUS: img/00000.png
img/00000.png: data/infected_00000.dat plot.py
	-rm -rf img
	mkdir img
	-python3 plot.py

.PHONY: plot
plot: img/00000.png
	@echo done

.PHONY: clean
clean:
	-rm -f $(OBJS) $(ASMS) $(DEPS) $(HDRS:.h=.h.gch) $(EXEC) *.out
	-rm -rf img
	-rm -rf data
	-mkdir img
	-mkdir data
	@echo done

.PHONY: debug
debug: $(DEPS) $(EXEC)
	@echo done

.PHONY: asm
asm: $(DEPS) $(ASMS)
	@echo done

.PHONY: depend
depend: $(DEPS)
	@echo done

.PHONY: headers
headers: $(HDRS:.h=.h.gch)
	@echo done

.PHONY: dox
dox: Doxyfile
	doxygen Doxyfile

$(EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.s: %.c
	$(CC) -S -g $(CFLAGS) -o $@ $<

%.d: %.c
	$(CC) $(DFLAGS) $< >$*.d

%.h.gch: %.h
	$(CC) -c $(CFLAGS) -o $@ $<

Doxyfile:
	doxygen -g
