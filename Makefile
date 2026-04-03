APP_NAME = hdbh

SRCDIR = src
SOURCES = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*.c)
OBJS = $(SOURCES:.cpp=.o)
OBJS := $(OBJS:.c=.o)

CXX = g++
CXXFLAGS = -Wall -O3 -std=c++20 -m64 -I. -Iinclude -Wno-unknown-pragmas
LDFLAGS = -lSDL2 -ldl -lpthread

DBGFLAGS = -Wall -O2 -g -std=c++20 -m64 -I. -Iinclude -Wno-unknown-pragmas

all: $(APP_NAME)

dbg: CXXFLAGS=$(DBGFLAGS)
dbg: clean $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(APP_NAME)