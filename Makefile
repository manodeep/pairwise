include common.mk

target:=pairwise_3d pairwise pairwise_3d_histogram

PAIRWISE_3D_SRC:=all_pairwise/pairwise_3d.c 
PAIRWISE_3D_OBJS:= $(PAIRWISE_3D_SRC:.c=.o)

ifeq ($(ISPC_AVAIL), 1)
  PAIRWISE_3D_ISPC_SRC:=all_pairwise/pairwise_3d.ispc
  PAIRWISE_3D_ISPC_OBJ:=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.o)
  PAIRWISE_3D_ISPC_HDR:=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.h)
endif

PAIRWISE_3D_HISTOGRAM_SRC:=pairwise_histogram/pairwise_3d_histogram.c
PAIRWISE_3D_HISTOGRAM_OBJ:=$(PAIRWISE_3D_HISTOGRAM_SRC:.c=.o)

PAIRWISE_SRC := all_pairwise/pairwise.c 
PAIRWISE_OBJ := $(PAIRWISE_SRC:.c=.o)

all: $(target) $(PAIRWISE_3D_SRC) $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile include/defs.h 

pairwise_3d:  $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_OBJS) include/defs.h utils.o progressbar.o 
	$(CC) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) utils.o progressbar.o -o $@ $(CLINK)

$(PAIRWISE_3D_ISPC_OBJ): $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile include/defs.h 
	ispc --addressing=64 $(OPT) -O3 $(PAIRWISE_3D_ISPC_SRC)  -h $(PAIRWISE_3D_ISPC_HDR) --target=avx1-i64x4 -o $@

$(PAIRWISE_3D_OBJS): $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $(PAIRWISE_3D_SRC) -o $@

pairwise: $(PAIRWISE_OBJ) utils.o common.mk Makefile include/defs.h progressbar.o 
	$(CC) $(PAIRWISE_OBJ) utils.o progressbar.o -o $@ $(CLINK) $(BLAS_LINK)

$(PAIRWISE_OBJ): $(PAIRWISE_SRC) common.mk Makefile include/defs.h
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) $(BLAS_INCLUDE) -c $< -o $@ 

$(PAIRWISE_3D_HISTOGRAM_OBJ): $(PAIRWISE_3D_HISTOGRAM_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $< -o $@

pairwise_3d_histogram: $(PAIRWISE_3D_HISTOGRAM_OBJ) utils.o include/defs.h progressbar.o 
	$(CC) utils.o progressbar.o $< -o $@ $(CLINK)

utils.o: utils/utils.c utils/utils.h
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

progressbar.o: utils/progressbar.c utils/progressbar.h
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

.PHONY: clean clena celan install lib tests distclean

clean:
	$(RM) $(target) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_ISPC_HDR) $(PAIRWISE_3D_HISTOGRAM_OBJ) pairwise_3d_histogram utils.o progressbar.o 

clena: clean

celan: celan
