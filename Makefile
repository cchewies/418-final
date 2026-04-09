APP_STD = hdbh
APP_NR  = hdbh_norender
APP_ECE = hdbh_ece

SRCDIR = src
SOURCES = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*.c)

# Build dirs
BUILD_STD = build/standard
BUILD_NR  = build/norender
BUILD_ECE = build/ece

# Compilers
MPICXX = mpic++
GXX = g++

# Flags
COMMON_FLAGS = -Wall -O3 -std=c++2a -m64 -I. -Iinclude -Wno-unknown-pragmas
FULL_LIBS = -lSDL2 -ldl -lpthread
NO_SDL_LIBS = -ldl -lpthread

# Object lists
OBJS_STD = $(patsubst $(SRCDIR)/%, $(BUILD_STD)/%, $(SOURCES:.cpp=.o))
OBJS_STD := $(OBJS_STD:.c=.o)

OBJS_NR = $(patsubst $(SRCDIR)/%, $(BUILD_NR)/%, $(SOURCES:.cpp=.o))
OBJS_NR := $(OBJS_NR:.c=.o)

OBJS_ECE = $(patsubst $(SRCDIR)/%, $(BUILD_ECE)/%, $(SOURCES:.cpp=.o))
OBJS_ECE := $(OBJS_ECE:.c=.o)

all: standard norender ece

# --- STANDARD ---
standard: $(APP_STD)

$(APP_STD): $(OBJS_STD)
	$(MPICXX) $(COMMON_FLAGS) -DRENDER_ENABLED -DUSE_MPI -o $@ $^ $(FULL_LIBS)

# --- NO RENDER ---
norender: $(APP_NR)

$(APP_NR): $(OBJS_NR)
	$(MPICXX) $(COMMON_FLAGS) -DUSE_MPI -o $@ $^ $(FULL_LIBS)

# --- ECE ---
ece: $(APP_ECE)

$(APP_ECE): $(OBJS_ECE)
	$(GXX) $(COMMON_FLAGS) -o $@ $^ $(NO_SDL_LIBS)

# --- Compile rules ---
$(BUILD_STD)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(MPICXX) $(COMMON_FLAGS) -DRENDER_ENABLED -DUSE_MPI -c $< -o $@

$(BUILD_STD)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(MPICXX) $(COMMON_FLAGS) -DRENDER_ENABLED -DUSE_MPI -c $< -o $@

$(BUILD_NR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(MPICXX) $(COMMON_FLAGS) -DUSE_MPI -c $< -o $@

$(BUILD_NR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(MPICXX) $(COMMON_FLAGS) -DUSE_MPI -c $< -o $@

$(BUILD_ECE)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(GXX) $(COMMON_FLAGS) -c $< -o $@

$(BUILD_ECE)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(GXX) $(COMMON_FLAGS) -c $< -o $@

clean:
	rm -rf build $(APP_STD) $(APP_NR) $(APP_ECE)