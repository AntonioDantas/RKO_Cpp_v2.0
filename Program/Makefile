# Detect the operating system
OS_RAW := $(shell uname -s)

# Verifica se a string retornada contém "MINGW", "CYGWIN" ou "MSYS"
# Isso cobre quase todos os ambientes de compilação no Windows
ifneq (,$(findstring MINGW,$(OS_RAW)))
    OS := Windows_NT
else ifneq (,$(findstring CYGWIN,$(OS_RAW)))
    OS := Windows_NT
else ifneq (,$(findstring MSYS,$(OS_RAW)))
    OS := Windows_NT
else
    OS := $(OS_RAW)
endif

# Configurações baseadas no SO detectado
ifeq ($(OS), Windows_NT)
    # No ambiente MinGW/MSYS, o 'rm' geralmente funciona melhor que o 'del'
    # Mas mantive compatível com o seu original, ajustando se necessário.
    RM := rm -f
    EXE := .exe
    CCC := g++
else ifeq ($(OS), Linux)
    RM := rm -f
    EXE :=
    CCC := g++
else ifeq ($(OS), Darwin)  
    RM := rm -f
    EXE :=
    CCC := g++-14
else
    # Fallback/Padrão caso a detecção falhe (Evita o erro "CCC vazio")
    RM := rm -f
    EXE :=
    CCC := g++
endif

TARGETDIR = .
CCFLAGS = -fopenmp -g -O3 -Wall -std=c++2a -Wno-unused-but-set-variable
TARGET = $(TARGETDIR)/runTest$(EXE)

debug: CCFLAGS += -g
release: CCFLAGS += -march=native

OBJS = $(TARGETDIR)/src/Main/main.o

# Rule to generate the final binary
$(TARGET): $(OBJS)
	$(CCC) $(CCFLAGS) -o $@ $(OBJS)

# Compilation of each object file
$(TARGETDIR)/src/Main/main.o: $(TARGETDIR)/src/Main/main.cpp
	$(CCC) $(CCFLAGS) -c $< -o $@

# Rule for running tests
test: $(TARGET)
	./$(TARGET)

# Cleaning rule
clean:
	$(RM) $(OBJS) $(TARGET)

# Rule for compiling in debug mode
debug: clean
	$(MAKE) $(TARGET) CCC=$(CCC)

# Rule to recompile
rebuild: 
	$(MAKE) debug
