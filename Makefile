CC=gcc

CFLAGS=-Wall -O2

OBJSti= 	ti_distance.o print_matrix.o make_rotation_matrix.o \
		rotate_tensor.o norm_matrix.o matrix_times_vector.o \
		vector_to_angles.o read_matrix.o find_ti.o

OBJSortho= 	print_matrix.o make_rotation_matrix.o rotate_tensor.o \
		norm_matrix.o matrix_times_vector.o \
		vector_to_angles.o ortho_distance.o quaternion_to_matrix.o \
		ti_distance.o read_matrix.o find_ortho.o

all: titest orthotest

clean:
	\rm titest orthotest *.o

titest: $(OBJSti) titest.o 
	gcc $(CFLAGS) titest.o $(OBJSti) -o $@ -lm -static

orthotest: $(OBJSortho) orthotest.o
	gcc $(CFLAGS) orthotest.o $(OBJSortho) -o $@ -lm -static
