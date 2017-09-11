# define compiler and path of libs
CPP = g++
CPPFLAG = -o
# define variable path
SRC = src
INC = include
BIN = bin
OBJ = obj
vpath %.cpp $(SRC)
vpath %.h $(INC)
vpath %.hpp $(INC)
vpath %.o $(OBJ)

all : bin/net.out bin/nets.out bin/spike.out bin/spike2.out bin/lfp.out bin/mi.out bin/mi_dd.out bin/mi_bd.out bin/ac.out bin/surrogate.out bin/surrogate2d.out bin/means.out bin/stds.out bin/transpose.out bin/split.out move
	echo 'All targets are done'

bin/net.out : get-config.o neuron.o	network.o connectivity_matrix.o main_net.o
	$(CPP) $(CPPFLAG) $@ $^

bin/nets.out : get-config.o neuron.o	network.o connectivity_matrix.o main_nets.o
	$(CPP) $(CPPFLAG) $@ $^

bin/spike.out : spike.o main_spike.o
	$(CPP) $(CPPFLAG) $@ $^

bin/spike2.out : spike.o main_spike2d.o
	$(CPP) $(CPPFLAG) $@ $^

bin/lfp.out : lfp.o main_lfp.o
	$(CPP) $(CPPFLAG) $@ $^

bin/mi.out : mi_uniform.o main_mi.o
	$(CPP) $(CPPFLAG) $@ $^

bin/mi_dd.out : mi_uniform.o main_mi_dd.o
	$(CPP) $(CPPFLAG) $@ $^

bin/mi_bd.out : mi_uniform.o main_mi_bd.o
	$(CPP) $(CPPFLAG) $@ $^

bin/ac.out : stationary.o main_ac.o
	$(CPP) $(CPPFLAG) $@ $^

bin/surrogate.out : surrogate.o main_surrogate.o
	$(CPP) $(CPPFLAG) $@ $^

bin/surrogate2d.out : surrogate.o main_surrogate2d.o
	$(CPP) $(CPPFLAG) $@ $^

bin/means.out : stationary.o main_means.o
	$(CPP) $(CPPFLAG) $@ $^

bin/stds.out : stationary.o main_stds.o
	$(CPP) $(CPPFLAG) $@ $^

bin/transpose.out : main_transpose.o
	$(CPP) $(CPPFLAG) $@ $^

bin/split.out : data_split.o
	$(CPP) $(CPPFLAG) $@ $^

.PHONY : move
move:
	mv *.o $(OBJ)
# .PHONY : clean
# clean:
# 	rm target $(OBJS)
