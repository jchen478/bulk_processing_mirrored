SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

COMPILER = icc 
EXE = stress_processing 

$(EXE): $(OBJS) 
	$(COMPILER) -o $@ $^

%.o: %.cpp
	$(COMPILER) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

