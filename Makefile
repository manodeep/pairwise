include common.mk

target=pairwise_3d pairwise pairwise_3d_histogram

PAIRWISE_3D_SRC=all_pairwise/pairwise_3d.c 
PAIRWISE_3D_OBJS= $(PAIRWISE_3D_SRC:.c=.o)

PAIRWISE_3D_ISPC_SRC=all_pairwise/pairwise_3d.ispc
PAIRWISE_3D_ISPC_OBJ=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.o)
PAIRWISE_3D_ISPC_HDR=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.h)

PAIRWISE_3D_HISTOGRAM_SRC=pairwise_histogram/pairwise_3d_histogram.c
PAIRWISE_3D_HISTOGRAM_OBJ=$(PAIRWISE_3D_HISTOGRAM_SRC:.c=.o)

PAIRWISE_SRC = all_pairwise/pairwise.c 
PAIRWISE_OBJ = $(PAIRWISE_SRC:.c=.o)

all: $(target) $(PAIRWISE_3D_SRC) $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile defs.h 

pairwise_3d:  $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_OBJS) defs.h  utils.o 
	$(CC) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) utils.o  -o $@ $(CLINK)

$(PAIRWISE_3D_ISPC_OBJ): $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile defs.h 
	ispc $(OPT) -O3 $(PAIRWISE_3D_ISPC_SRC)  -h $(PAIRWISE_3D_ISPC_HDR) --target=avx1-i64x4 -o $@

$(PAIRWISE_3D_OBJS): $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $(PAIRWISE_3D_SRC) -o $@

pairwise: $(PAIRWISE_OBJ) utils.o common.mk Makefile defs.h 
	$(CC) $(PAIRWISE_OBJ) utils.o  -o $@ $(CLINK) $(BLAS_LINK)

$(PAIRWISE_OBJ): $(PAIRWISE_SRC) common.mk Makefile defs.h
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) $(BLAS_INCLUDE) -c $< -o $@ 

$(PAIRWISE_3D_HISTOGRAM_OBJ): $(PAIRWISE_3D_HISTOGRAM_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $< -o $@

pairwise_3d_histogram: $(PAIRWISE_3D_HISTOGRAM_OBJ) utils.o defs.h 
	$(CC)  utils.o $< -o $@ $(CLINK)

utils.o: utils/utils.c utils/utils.h
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

.PHONY: clean clena celan install lib tests distclean

clean:
	$(RM) $(target) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_ISPC_HDR) 

clena: clean

celan: celan
