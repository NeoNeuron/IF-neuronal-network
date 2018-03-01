# define compiler and path of libs
CPP = g++
CPPFLAG = -o
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

all : net.out nets.out spike.out spike2.out lfp.out mi.out mi_dd.out mi_dd_LFP.out mi_bd.out mi_bb.out autocov.out surrogate.out surrogate2d.out means.out stds.out transpose.out split.out sta.out move
	echo 'All targets are done'

net.out : get-config.o neuron.o	network.o connectivity_matrix.o main_net.o
	$(CPP) $(CPPFLAG) $@ $^

nets.out : get-config.o neuron.o	network.o connectivity_matrix.o main_nets.o
	$(CPP) $(CPPFLAG) $@ $^

spike.out : spike.o main_spike.o
	$(CPP) $(CPPFLAG) $@ $^

spike2.out : spike.o main_spike2d.o
	$(CPP) $(CPPFLAG) $@ $^

lfp.out : lfp.o main_lfp.o
	$(CPP) $(CPPFLAG) $@ $^

mi.out : mi_uniform.o main_mi.o
	$(CPP) $(CPPFLAG) $@ $^

mi_dd.out : mi_uniform.o main_mi_dd.o
	$(CPP) $(CPPFLAG) $@ $^

mi_dd_LFP.out : mi_uniform.o main_mi_dd_LFP.o
	$(CPP) $(CPPFLAG) $@ $^

mi_bd.out : mi_uniform.o main_mi_bd.o
	$(CPP) $(CPPFLAG) $@ $^

mi_bb.out : mi_uniform.o main_mi_bb.o
	$(CPP) $(CPPFLAG) $@ $^

autocov.out : stationary.o main_autocov.o
	$(CPP) $(CPPFLAG) $@ $^

surrogate.out : surrogate.o main_surrogate.o
	$(CPP) $(CPPFLAG) $@ $^

surrogate2d.out : surrogate.o main_surrogate2d.o
	$(CPP) $(CPPFLAG) $@ $^

means.out : stationary.o main_means.o
	$(CPP) $(CPPFLAG) $@ $^

stds.out : stationary.o main_stds.o
	$(CPP) $(CPPFLAG) $@ $^

transpose.out : main_transpose.o
	$(CPP) $(CPPFLAG) $@ $^

split.out : data_split.o
	$(CPP) $(CPPFLAG) $@ $^

sta.out : sta.o
	$(CPP) $(CPPFLAG) $@ $^

.PHONY : move
move:
	mv *.o $(OBJ)
	mv *.out $(BIN)
# .PHONY : clean
# clean:
# 	rm target $(OBJS)
