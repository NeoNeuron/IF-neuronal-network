# define compiler and path of libs
CPPFLAGS = --std=c++11 -Wall -g -O2 
#LDLIBS = -lboost_program_options
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

.PHONY : all
all : net.out spike.out spike2.out lfp.out mi.out mi_dd.out mi_dd_LFP.out mi_bd.out mi_bd_2bins.out mi_bb.out mi_bd_unity.out transpose.out sta.out move
	echo 'All targets are done'

neuron.out : get-config.o neuron.o poisson_generator.o math_helper.o main_neuron.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

net.out : get-config.o neuron.o poisson_generator.o	network.o math_helper.o main_net.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

spike.out : spike.o main_spike.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

spike2.out : spike.o main_spike2d.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

lfp.out : get-config.o lfp.o main_lfp.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi.out : mi_uniform.o main_mi.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_dd.out : mi_uniform.o main_mi_dd.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_dd_LFP.out : mi_uniform.o main_mi_dd_LFP.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_bd.out : mi_uniform.o main_mi_bd.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_bd_2bins.out : mi_uniform.o main_mi_bd_2bins.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_bb.out : mi_uniform.o main_mi_bb.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

mi_bd_unity.out : spike.o lfp.o mi_uniform.o mi_bd_unity.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

transpose.out : main_transpose.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

sta.out : sta.o
	$(CXX) -MM $(CPPFLAGS) -o $@ $^

## Generate prerequisites
#%.d: %.cpp
#	@set -e; rm -f $@; \
#		$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
#		sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
#		rm -f $@.$$$$
#
#SRCDEP = $(SRCS:.cpp=.d)
#include $(SRCDEP)

.PHONY : move
move:
	mkdir -p $(OBJ)
	mv *.o $(OBJ)
	mkdir -p $(BIN)
	mv *.out $(BIN)

.PHONY : clean
clean:
	rm -rf $(OBJ)
