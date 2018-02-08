SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

COMPILER = gcc 
EXE = diffusion 
CFLAG = -lm

$(EXE): $(OBJS) 
	$(COMPILER) $(CFLAG) -o $@ $^

%.o: %.cpp
	$(COMPILER) $(CFLAG) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

