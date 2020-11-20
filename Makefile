
DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

COMPILER       = gcc
DEBUG          =
OPTIMIZATION   = -O3
FLOWTRACE      =
OPTIONS        = -std=c++11
LINK_OPTIONS   =
CFLAGS         = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)

LIBS     = -L /usr/local/lib -lm -lgsl -lgslcblas -lc++
INCLUDES = -I /usr/local/include -I rhic/include

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
OBJ = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE = cpu_vah

$(EXE): $(OBJ)
	echo "Linking: $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi


.SILENT:
