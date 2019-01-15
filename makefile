# define compiler and path of libs
CPP = g++
CPPFLAGS = --std=c++11 -g -O2 
# define variable path
SRC = src
INC = include
BIN = bin
OBJ = obj
BIN = bin
vpath %.cpp $(SRC)
vpath %.h $(INC)
vpath %.hpp $(INC)
vpath %.o $(OBJ)
vpath %.out $(BIN)

all : net.out spike.out spike2.out lfp.out mi.out mi_dd.out mi_dd_LFP.out mi_bd.out mi_bd_2bins.out mi_bb.out mi_bd_unity.out transpose.out sta.out move
	echo 'All targets are done'

neuron.out : get-config.o neuron.o poisson_generator.o math_helper.o main_neuron.o
	$(CPP) $(CPPFLAGS) -o $@ $^

net.out : get-config.o neuron.o poisson_generator.o	network.o math_helper.o main_net.o
	$(CPP) $(CPPFLAGS) -o $@ $^

spike.out : spike.o main_spike.o
	$(CPP) $(CPPFLAGS) -o $@ $^

spike2.out : spike.o main_spike2d.o
	$(CPP) $(CPPFLAGS) -o $@ $^

lfp.out : get-config.o lfp.o main_lfp.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi.out : mi_uniform.o main_mi.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_dd.out : mi_uniform.o main_mi_dd.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_dd_LFP.out : mi_uniform.o main_mi_dd_LFP.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_bd.out : mi_uniform.o main_mi_bd.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_bd_2bins.out : mi_uniform.o main_mi_bd_2bins.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_bb.out : mi_uniform.o main_mi_bb.o
	$(CPP) $(CPPFLAGS) -o $@ $^

mi_bd_unity.out : spike.o lfp.o mi_uniform.o mi_bd_unity.o
	$(CPP) $(CPPFLAGS) -o $@ $^

transpose.out : main_transpose.o
	$(CPP) $(CPPFLAGS) -o $@ $^

sta.out : sta.o
	$(CPP) $(CPPFLAGS) -o $@ $^

.PHONY : move
move:
	mkdir -p $(OBJ)
	mv *.o $(OBJ)
	mkdir -p $(BIN)
	mv *.out $(BIN)

.PHONY : clean
clean:
	rm -rf $(OBJ)
