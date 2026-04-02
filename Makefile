APP_NAME = hdbh

SRCDIR = src
OBJS = $(SRCDIR)/entry.o

CXX = mpic++
CXXFLAGS = -Wall -O3 -std=c++20 -m64 -I. -fopenmp -Wno-unknown-pragmas
DBGFLAGS = -Wall -O2 -g -std=c++20 -m64 -I. -fopenmp -Wno-unknown-pragmas

all: $(APP_NAME)
dbg: CXXFLAGS=$(DBGFLAGS)
dbg: clean $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(APP_NAME)