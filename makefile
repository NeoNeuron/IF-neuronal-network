# define compiler and path of libs
CPPFLAGS = --std=c++11 -Wall
CXXFLAGS = -g -O2
#LDLIBS = -lboost_program_options
# define variable path
DIR_SRC = src
DIR_INC = include
DIR_BIN = bin
DIR_OBJ = obj
DIR_DEP = dep
DIRS = $(DIR_BIN) $(DIR_DEP) $(DIR_OBJ)
vpath %.cpp $(DIR_SRC)
vpath %.h 	$(DIR_INC)
vpath %.hpp $(DIR_INC)
vpath %.d 	$(DIR_DEP)
vpath %.o 	$(DIR_OBJ)
vpath %.out $(DIR_BIN)
# include all source in DIR_SRC
SRCS := $(notdir $(wildcard $(DIR_SRC)/*.cpp))
DEPS = $(SRCS:.cpp=.d)
DEPS := $(addprefix $(DIR_DEP)/, $(DEPS))
#OBJS := $(SRCS:.cpp=.o)
#OBJS := $(addprefix $(DIR_OBJ)/, $(OBJS))

.PHONY : all
all : $(DIRS) net.out spike.out lfp.out mi.out mi_bd.out mi_bd_2bins.out move
	@echo 'All targets are done'

neuron.out : get-config.o neuron.o poisson_generator.o math_helper.o main_neuron.o
	$(CXX) -o $@ $^

net.out : get-config.o neuron.o poisson_generator.o	network.o math_helper.o main_net.o
	$(CXX) -o $@ $^

spike.out : spike.o main_spike.o
	$(CXX) -o $@ $^

spike2.out : spike.o main_spike2d.o
	$(CXX) -o $@ $^

lfp.out : get-config.o lfp.o main_lfp.o
	$(CXX) -o $@ $^

mi.out : mi_uniform.o main_mi.o
	$(CXX) -o $@ $^

mi_dd.out : mi_uniform.o main_mi_dd.o
	$(CXX) -o $@ $^

mi_dd_LFP.out : mi_uniform.o main_mi_dd_LFP.o
	$(CXX) -o $@ $^

mi_bd.out : mi_uniform.o main_mi_bd.o
	$(CXX) -o $@ $^

mi_bd_2bins.out : mi_uniform.o main_mi_bd_2bins.o
	$(CXX) -o $@ $^

mi_bb.out : mi_uniform.o main_mi_bb.o
	$(CXX) -o $@ $^

mi_bd_unity.out : spike.o lfp.o mi_uniform.o mi_bd_unity.o
	$(CXX) -o $@ $^

transpose.out : main_transpose.o
	$(CXX) -o $@ $^

sta.out : sta.o
	$(CXX) -o $@ $^

ifeq ($(wildcard $(DIR_DEP)),)
DEP_DIR_DEP := $(DIR_DEP)
endif

$(DIRS) : 
	mkdir -p $@

# Generate prerequisites
$(DIR_DEP)/%.d : $(DEP_DIR_DEP) %.cpp
	@echo "Making $@ ..."
	@set -e; rm -f $@; \
		$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
		sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
		rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

#$(DIR_OBJ)/%.o : $(DIR_OBJ) %.cpp
#	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $^

.PHONY : move
move:
	@mv *.o $(DIR_OBJ)
	@mv *.out $(DIR_BIN)

.PHONY : clean
clean:
	rm -rf $(DIR_OBJ) $(DIR_BIN) $(DIR_DEP)
