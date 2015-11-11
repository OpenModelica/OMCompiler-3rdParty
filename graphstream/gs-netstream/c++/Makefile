BINDIR=bin
OBJDIR=obj
SRCDIR=src
INCDIR=include

ifeq (darwin,$(findstring darwin,$(shell $(CC) -v 2>&1 | grep ^Target)))
	AR=libtool
	ARFLAGS=-static -o
else
	AR=ar
	ARFLAGS=rcs
endif

CWARN = -W -Wall -Wshadow -Wimplicit -Wreturn-type -Wcomment -Wtrigraphs -Wformat -Wparentheses -Wpointer-arith -Wuninitialized -O

CDBG = -g $(CWARN) -fno-inline
DFLAGS = -I$(INCDIR) -g $(CWARN) -fno-inline -DDEBUG=1

ifeq (mingw,$(findstring mingw,$(shell $(CC) -v 2>&1 | grep ^Target)))
	SOCKET_LIB = -lwsock32
else
	SOCKET_LIB =
endif


CTAG = ctags
CTAGFILE = filelist
# src, object and bin files
OUT_BINARY := netstream-main
OUT_LIBRARY := libnetstream.a

HEADERS = $(INCDIR)/global.h

OBJS = \
    $(OBJDIR)/netstream-storage.o  \
	$(OBJDIR)/netstream-socket.o \
    $(OBJDIR)/netstream-sender.o \
	$(OBJDIR)/netstream-main.o

.SECONDARY:

#-- Rules
all: $(OUT_BINARY) $(OUT_LIBRARY)

$(OUT_BINARY): $(OBJS)
	$(CXX) $(CFLAGS) -I$(INCDIR) $(OBJS) -o $(OUT_BINARY) $(SOCKET_LIB)

$(OUT_LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_LIBRARY)  $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c $? -o $@

.PHONY: clean depend fresh

tag:
	find src/*.cpp include/*.h > filelist
	$(CTAG) -L $(CTAGFILE)

clean:
	-rm -f $(OBJDIR)/*.o $(PARSE_C) $(PARSE_H)
	-rm -f $(OBJDIR)/*/*.o $(PARSE_C) $(PARSE_H)
	-rm -f $(SRCDIR)/*.output $(LEX_C)
	-rm -f */*~ *~ core
	-rm -f $(BINDIR)/*
	-rm -f $(OUT_BINARY) $(OUT_BINARY).exe
	-rm -f $(OUT_LIBRARY)
	-rm -rf $(OBJDIR)

fresh: clean all
