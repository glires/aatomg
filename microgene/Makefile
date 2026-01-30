CC = gcc
CFLAGS = -O3 -std=c99 -Wall
BINDIR = ./bin
SRCDIR = ./src

TARGETS = $(BINDIR)/gor $(BINDIR)/aatobase $(BINDIR)/edmg0 $(BINDIR)/edmg1 \
          $(BINDIR)/edmg2 $(BINDIR)/microgene $(BINDIR)/ssscore \
          $(BINDIR)/uniqaa $(BINDIR)/scoreab

all: $(TARGETS)

$(BINDIR)/gor: $(SRCDIR)/gor.c
	$(CC) $(CFLAGS) -o $@ $< -lm

$(BINDIR)/aatobase: $(SRCDIR)/aatobase.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/edmg0: $(SRCDIR)/edmg0.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/edmg1: $(SRCDIR)/edmg1.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/edmg2: $(SRCDIR)/edmg2.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/microgene: $(SRCDIR)/microgene.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/ssscore: $(SRCDIR)/ssscore.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/uniqaa: $(SRCDIR)/uniqaa.c
	$(CC) $(CFLAGS) -o $@ $<

$(BINDIR)/scoreab: $(SRCDIR)/scoreab.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)

.PHONY: all clean
