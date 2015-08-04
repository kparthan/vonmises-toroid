#******************************************************************************
# TARGETS : ./bin/phipsiomega
# AUTHOR  : ARUN S KONAGURTHU
#******************************************************************************
#directories
phipsiomega = .

#local directories
SRC = $(phipsiomega)/src
OBJ = $(phipsiomega)/obj
BIN = $(phipsiomega)/bin

#compiler options
CPP = g++ 
CPPFLAGS = -O3 
LDFLAGS =
.cpp.o:
	$(CPP) $(CPPFLAGS) -c -o $@ $<

#macros
OBJECTS = $(OBJ)/pdbClass.o $(OBJ)/superpose3D.o  $(OBJ)/geometry3D.o $(OBJ)/main.o
ALL = $(BIN)/phipsiomega

#targets
all: $(ALL)


#------------------------------------------------------------------------------
$(BIN)/phipsiomega: $(OBJECTS) 
	$(CPP) $(CPPFLAGS)  $(LDFLAGS) -o $@ $(OBJECTS)

$(OBJ)/pdbClass.o: $(SRC)/pdbClass.h $(SRC)/pdbClass.cpp $(SRC)/geometry3D.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/pdbClass.cpp -o $(OBJ)/pdbClass.o
		
$(OBJ)/superpose3D.o: $(SRC)/superpose3D.h $(SRC)/superpose3D.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/superpose3D.cpp -o $(OBJ)/superpose3D.o

$(OBJ)/geometry3D.o: $(SRC)/geometry3D.h $(SRC)/geometry3D.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/geometry3D.cpp -o $(OBJ)/geometry3D.o
		
$(OBJ)/main.o: $(SRC)/main.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/main.cpp -o $(OBJ)/main.o

#------------------------------------------------------------------------------
clean:
	rm $(OBJ)/*.o

remove: clean
	rm  $(BIN)/*
#------------------------------------------------------------------------------

#******************************************************************************
