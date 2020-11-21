DIR_MAIN = ./
DIR_SRC = ./
DIR_OBJ = ./

COMPILER = icpc
OPTIMIZATION = -O3

CFLAGS= $(OPTIMIZATION) $(OPTIONS)
CFLAGS = $(CFLAG) -Wno-comment

LIBS = -lm

INCLUDES =

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
OBJ =$(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE = data

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

clean:
	echo "\nDeleting executable and object files\n"
	rm -rf $(EXE)
	rm *.o

.SILENT:

