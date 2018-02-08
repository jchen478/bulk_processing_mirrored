SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

COMPILER = gcc -lm 
EXE = elastic_processing 

$(EXE): $(OBJS) 
	$(COMPILER) -o $@ $^

%.o: %.cpp
	$(COMPILER) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

