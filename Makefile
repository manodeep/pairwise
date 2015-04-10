include common.mk

target=pairwise_3d pairwise pairwise_3d_histogram

PAIRWISE_3D_SRC=all_pairwise/pairwise_3d.c 
PAIRWISE_3D_OBJS= $(PAIRWISE_3D_SRC:.c=.o)

PAIRWISE_3D_ISPC_SRC=all_pairwise/pairwise_3d.ispc
PAIRWISE_3D_ISPC_OBJ=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.o)
PAIRWISE_3D_ISPC_HDR=$(PAIRWISE_3D_ISPC_SRC:.ispc=_ispc.h)

PAIRWISE_3D_HISTOGRAM_SRC=pairwise_histogram/pairwise_3d_histogram.c
PAIRWISE_3D_HISTOGRAM_OBJ=$(PAIRWISE_3D_HISTOGRAM_SRC:.c=.o)

all: $(target) $(PAIRWISE_3D_SRC) $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile

pairwise_3d:  $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_OBJS)
	$(CC) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) $(CLINK) -o $@

$(PAIRWISE_3D_ISPC_OBJ): $(PAIRWISE_3D_ISPC_SRC) common.mk Makefile
	ispc $(OPT) -O3 $(PAIRWISE_3D_ISPC_SRC)  -h $(PAIRWISE_3D_ISPC_HDR) --target=avx1-i64x4 -o $@

$(PAIRWISE_3D_OBJS): $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $(PAIRWISE_3D_SRC) -o $@

pairwise: all_pairwise/pairwise.c common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) $(BLAS_INCLUDE)  $< $(CLINK) $(BLAS_LINK) -o $@


$(PAIRWISE_3D_HISTOGRAM_OBJ): $(PAIRWISE_3D_HISTOGRAM_SRC) common.mk Makefile
	$(CC) $(OPT) $(INCLUDE) $(CFLAGS) -c $< -o $@

pairwise_3d_histogram: $(PAIRWISE_3D_HISTOGRAM_OBJ) utils.o
	$(CC) $(CLINK) utils.o $< -o $@

utils.o: utils/utils.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

.PHONY: clean clena celan install lib tests distclean

clean:
	$(RM) $(target) $(PAIRWISE_3D_OBJS) $(PAIRWISE_3D_ISPC_OBJ) $(PAIRWISE_3D_ISPC_HDR) 

clena: clean

celan: celan
