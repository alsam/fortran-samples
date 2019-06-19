.SUFFIXES: .c .o .h .cuf .f90

#F90FLAGS = -fast
F90FLAGS = -g  # -fimplicit-none
CUDAFLAGS = -Mcuda=cc20,cuda3.2,keepptx
CUDAFLAGS = -rc=rc4.0 -Mcuda=cc20,keepptx
 
%.o: %.f90
	gfortran $(F90FLAGS) -c -o $@ $<

%.o: %.cuf
	pgf90 $(CUDAFLAGS) $(F90FLAGS) -c -o $@ $<

%.o: %.c
	pgcc -c -o $@ $<

# Select 32- or 64-bit libs
#CUDALIB=/usr/local/cuda/lib
#CUDALIB=/usr/local/cuda/lib64
CUDALIB=/usr/local/cuda-tot/lib64

#TARGET = helm2d
TARGET = body_ax

OBJS = 	precision_mod.o parameters_mod.o common_mod01.o common_mod02.o gauss_legendre.o \
	body_ax_geo.o  	body_ax_sdlp.o  elliptic_integral.o  elm_line.o lgf_ax_fs.o \
	body_ax.o      	body_ax_vel.o   elm_arc.o  \
	gel.o        	lvr_fs.o

all: $(TARGET)

#$(TARGET):  $(OBJS)
#	pgf90 $(CUDAFLAGS) $(F90FLAGS) -o $@ $^ -L$(CUDALIB) -lcufft -lcublas

$(TARGET): $(OBJS)
	gfortran $(F90FLAGS) -o $@ $^


utest:
	echo '51' | ./$(TARGET)
	diff linear_system.solution.text linear_system.solution.text.golden
	diff body_ax.str body_ax.str.golden

#
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJS) body_ax *.mod
	rm -f body_ax.out body_ax.str xyunit
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJS) *.mod
	rm -f body_ax.out body_ax.str xyunit
