SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

COMPILER = g++ 
CFLAGS = -O2
EXE = Slice 

$(EXE): $(OBJS) 
	$(COMPILER) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(COMPILER) $(CFLAGS) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

