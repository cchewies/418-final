APP_NAME = hdbh

SRCDIR = src
SOURCES = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*.c)
OBJS = $(SOURCES:.cpp=.o)
OBJS := $(OBJS:.c=.o)

CXX = mpic++
CXXFLAGS = -Wall -O3 -std=c++20 -m64 -I. -Iinclude -Wno-unknown-pragmas
LDFLAGS = -lSDL2 -ldl -lpthread

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(APP_NAME)