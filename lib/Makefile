#
# Config
#

CMACH = PC_LINUX1

LIB_DIR = $(HOME)/etc/rhdf5/lib
MOD_DIR = $(HOME)/etc/rhdf5/mod
ARC = librhdf5_utils.a
MOD = rhdf5_utils.mod
DEST_ARC = $(LIB_DIR)/$(ARC)
DEST_MOD = $(MOD_DIR)/$(MOD)

INCLUDES = -I/Users/steve/hdf5-1.8.6-macosx64-static/include

OBJS = $(ARC)(rhdf5_f2c.o) \
       $(ARC)(rhdf5_utils.o)

ARC = librhdf5_utils.a

F_COMP = pgf90
F_OPTS = -Munroll -Mnoframe -O2 -pc 64 -tp x64 -Mfree -fastsse -time
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

C_COMP = gcc
C_OPTS = -O3 -DUNDERSCORE -DLITTLE -D$(CMACH)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

ARCHIVE = ar rs

#
# Targets
#

all: $(DEST_ARC) $(DEST_MOD)

$(DEST_ARC): $(ARC)
	cp $(ARC) $(DEST_ARC)

$(DEST_MOD): $(MOD)
	cp $(MOD) $(DEST_MOD)

$(ARC): $(OBJS)

clean:
	rm -f *.a
	rm -f *.o
	rm -f *.mod

#
# Rules
#

.SUFFIXES: .f90 .c .o .a

# f90 rule.

.f90.a:
	@echo ""
	$(F_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.f90=.o)
	rm -f $(<F:.f90=.o)

# c rule.

.c.a:
	@echo ""
	$(C_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.c=.o)
	rm -f $(<F:.c=.o)

