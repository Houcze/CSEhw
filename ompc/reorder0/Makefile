CC     = mpicc
OPTFLAGS= -o3 -mprefer-vector-width=512
TARGET    = L3-rsol0.exe

.SUFFIXES:
.SUFFIXES: .o .c 

.c.o:
	$(CC) -g -c $(CFLAGS) $(OPTFLAGS) $< -o $@ -lm -lgomp

OBJS = \
	input.o \
	pointer_init.o \
	boundary_cell.o \
	cell_metrics.o \
	mc.o \
	cm.o \
	rcm.o \
	cmrcm.o \
	poi_gen.o \
	solver_ICCG_mc.o \
	solver_ICCG_mc_ft.o \
	outucd.o \
	allocate.o \
	main.o

HEADERS = \
	struct.h \
	struct_ext.h \
	pcg.h \
	pcg_ext.h \
	input.h \
	pointer_init.h \
	boundary_cell.h \
	cell_metrics.h \
	mc.h \
	cm.h \
	rcm.h \
	cmrcm.h \
	poi_gen.h \
	solver_ICCG_mc.h \
	solver_ICCG_mc_ft.h \
	allocate.h \
	outucd.h

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) $(OPTFLAGS) -o $@ $(OBJS) $(LDFLAGS) -lm -lgomp

$(OBJS): $(HEADERS)

clean:
	rm -f *.o $(TARGET) *.log *~
