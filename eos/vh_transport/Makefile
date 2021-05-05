
UNAME := $(shell uname)
NVCC := $(shell command -v nvcc 2> /dev/null)

DIR_MAIN = ./
DIR_SRC = ./
DIR_OBJ = ./

CFLAGS= $(OPTIMIZATION) $(OPTIONS)


# Compiler options
ifdef NVCC
COMPILER = nvcc
OPTIMIZATION := -O3
OPTIONS := -Wno-deprecated-gpu-targets
endif
ifndef NVCC
COMPILER = gcc
OPTIMIZATION = -O3
CFLAGS := $(CFLAG) -Wno-comment
endif

# Library options
ifeq ($(UNAME), Linux)
LIBS := -lm -lgsl
endif
ifeq ($(UNAME), Darwin)
LIBS := -lc++ -lm
#LIBS := -lm -lgsl
endif

INCLUDES = -I include

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
OBJ =$(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE = vh_transport

# $@ = target
# $< = dependency
# $^ = + everything else

$(EXE): $(OBJ)
	echo "\nLinking   $@ ($(COMPILER))"
	$(COMPILER) -o $@ $^ $(LIBS) $(INCLUDES)
	echo "\nRunning...\n"
	$(DIR_MAIN)$(EXE)
	rm *.o

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	echo "Compiling $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

run:
	echo "\nRunning...\n"
	$(DIR_MAIN)$(EXE)

clean:
	echo "\nDeleting executable and object files\n"
	rm -rf $(EXE)
	rm *.o

.SILENT:

