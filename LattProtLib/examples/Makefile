#-----------------------------------------------------------------------
# Makefile for example programs of the lattice_prot_lib
# Uses lattice_prot_lib as source
# Generates config.hpp file from params.json file and then compiles it
# along with the rest of the source files to have compile time parameters
#-----------------------------------------------------------------------

# === User-defined paths ===
LIB_DIR        := /home/nicolas/quantevol\ Dropbox/statbio\ CIRB/Nicolas/codes/lattice_prot_lib
SCRIPT_DIR     := $(LIB_DIR)/pyScripts

# Config files
CONFIG_JSON_MH_PULL_CUBIC := config_files/params_MH_PULL_cubic.json
CONFIG_MH_PULL_CUBIC_HPP  := src/config_MH_PULL_cubic.hpp

CONFIG_JSON_MH_PULL_SQUARE := config_files/params_MH_PULL_square.json
CONFIG_MH_PULL_SQUARE_HPP  := src/config_MH_PULL_square.hpp

CONFIG_JSON_MH_VSHD_SQUARE := config_files/params_MH_VSHD_square.json
CONFIG_MH_VSHD_SQUARE_HPP  := src/config_MH_VSHD_square.hpp

CONFIG_JSON_SA_TM_SQUARE := config_files/params_SA_TM_square.json
CONFIG_SA_TM_SQUARE_HPP  := src/config_SA_TM_square.hpp

# === Compiler and flags ===
CXX            := g++
CXXFLAGS       := -O3 -march=native -funroll-loops -std=c++23 -Wall -Wextra \
                  -I$(LIB_DIR)/include \
                  -I/home/nicolas/HighFive/include
LDFLAGS        := -fopenmp $(HDF5_LIBS)

# === Detect HDF5 variant ===
ifneq (, $(shell ldconfig -p | grep libhdf5_serial_cpp))
    HDF5_LIBS    = -lhdf5_serial_cpp -lhdf5_serial
    HDF5_INCLUDE = -I/usr/include/hdf5/serial
else
    HDF5_LIBS    = -lhdf5_cpp -lhdf5
    HDF5_INCLUDE =
endif
CXXFLAGS += $(HDF5_INCLUDE)
LDFLAGS  += $(HDF5_LIBS)

# === Directories ===
SRC_DIR         := $(LIB_DIR)/src
PROJECT_SRC_DIR := src
OBJ_DIR         := obj
BIN_DIR         := bin

# === Library sources ===
LIB_SOURCES     := $(wildcard $(SRC_DIR)/*.cpp)
LIB_OBJECTS     := $(LIB_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# === Project sources (will filter later) ===
PROJECT_SOURCES := $(wildcard $(PROJECT_SRC_DIR)/*.cpp)


# ------------------------------------------------------------
# Targets
# ------------------------------------------------------------
.PHONY: all evo sim clean

all: MH_VSHD_square MH_PULL_square MH_PULL_cubic SA_TM_square WL_PULL_square

# ------------------------------------------------------------
# Build MH_PULL_square
# ------------------------------------------------------------
MH_PULL_square: $(CONFIG_MH_PULL_SQUARE_HPP) $(OBJ_DIR)/MH_PULL_square.o $(LIB_OBJECTS)
	@echo "[Link] MH_PULL_square"
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LIB_OBJECTS) $(OBJ_DIR)/MH_PULL_square.o -o $(BIN_DIR)/MH_PULL_square.out $(LDFLAGS)

$(CONFIG_MH_PULL_SQUARE_HPP): $(CONFIG_JSON_MH_PULL_SQUARE) $(SCRIPT_DIR)/generate_config.py
	@echo "[Make] Generating $@ from $(CONFIG_JSON_MH_PULL_SQUARE)"
	@mkdir -p $(dir $@)
	python3 $(SCRIPT_DIR)/generate_config.py $(CONFIG_JSON_MH_PULL_SQUARE) $(CONFIG_MH_PULL_SQUARE_HPP)

$(OBJ_DIR)/MH_PULL_square.o: $(PROJECT_SRC_DIR)/MH_PULL_square.cpp $(CONFIG_MH_PULL_SQUARE_HPP)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Build MH_PULL_cubic
# ------------------------------------------------------------
MH_PULL_cubic: $(CONFIG_MH_PULL_CUBIC_HPP) $(OBJ_DIR)/MH_PULL_cubic.o $(LIB_OBJECTS)
	@echo "[Link] MH_PULL_cubic"
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LIB_OBJECTS) $(OBJ_DIR)/MH_PULL_cubic.o -o $(BIN_DIR)/MH_PULL_cubic.out $(LDFLAGS)

$(CONFIG_MH_PULL_CUBIC_HPP): $(CONFIG_JSON_MH_PULL_CUBIC) $(SCRIPT_DIR)/generate_config.py
	@echo "[Make] Generating $@ from $(CONFIG_JSON_MH_PULL_CUBIC)"
	@mkdir -p $(dir $@)
	python3 $(SCRIPT_DIR)/generate_config.py $(CONFIG_JSON_MH_PULL_CUBIC) $(CONFIG_MH_PULL_CUBIC_HPP)

$(OBJ_DIR)/MH_PULL_cubic.o: $(PROJECT_SRC_DIR)/MH_PULL_cubic.cpp $(CONFIG_MH_PULL_CUBIC_HPP)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Build MH_VSHD_square
# ------------------------------------------------------------
MH_VSHD_square: $(CONFIG_MH_VSHD_SQUARE_HPP) $(OBJ_DIR)/MH_VSHD_square.o $(LIB_OBJECTS)
	@echo "[Link] MH_VSHD_square"
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LIB_OBJECTS) $(OBJ_DIR)/MH_VSHD_square.o -o $(BIN_DIR)/MH_VSHD_square.out $(LDFLAGS)

$(CONFIG_MH_VSHD_SQUARE_HPP): $(CONFIG_JSON_MH_VSHD_SQUARE) $(SCRIPT_DIR)/generate_config.py
	@echo "[Make] Generating $@ from $(CONFIG_JSON_MH_VSHD_SQUARE)"
	@mkdir -p $(dir $@)
	python3 $(SCRIPT_DIR)/generate_config.py $(CONFIG_JSON_MH_VSHD_SQUARE) $(CONFIG_MH_VSHD_SQUARE_HPP)

$(OBJ_DIR)/MH_VSHD_square.o: $(PROJECT_SRC_DIR)/MH_VSHD_square.cpp $(CONFIG_MH_VSHD_SQUARE_HPP)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Build WL_PULL_square
# ------------------------------------------------------------
WL_PULL_square: $(OBJ_DIR)/WL_PULL_square.o $(LIB_OBJECTS)
	@echo "[Link] WL_PULL_square"
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LIB_OBJECTS) $(OBJ_DIR)/WL_PULL_square.o -o $(BIN_DIR)/WL_PULL_square.out $(LDFLAGS)

$(OBJ_DIR)/WL_PULL_square.o: $(PROJECT_SRC_DIR)/WL_PULL_square.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Build SA_TM_square
# ------------------------------------------------------------
SA_TM_square: $(CONFIG_SA_TM_SQUARE_HPP) $(OBJ_DIR)/SA_TM_square.o $(LIB_OBJECTS)
	@echo "[Link] SA_TM_square"
	@mkdir -p $(BIN_DIR)
	$(CXX) $(LIB_OBJECTS) $(OBJ_DIR)/SA_TM_square.o -o $(BIN_DIR)/SA_TM_square.out $(LDFLAGS)

$(CONFIG_SA_TM_SQUARE_HPP): $(CONFIG_JSON_SA_TM_SQUARE) $(SCRIPT_DIR)/generate_config.py
	@echo "[Make] Generating $@ from $(CONFIG_JSON_SA_TM_SQUARE)"
	@mkdir -p $(dir $@)
	python3 $(SCRIPT_DIR)/generate_config.py $(CONFIG_JSON_SA_TM_SQUARE) $(CONFIG_SA_TM_SQUARE_HPP)

$(OBJ_DIR)/SA_TM_square.o: $(PROJECT_SRC_DIR)/SA_TM_square.cpp $(CONFIG_SA_TM_SQUARE_HPP)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Compile library sources
# ------------------------------------------------------------
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ------------------------------------------------------------
# Clean
# ------------------------------------------------------------
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(CONFIG_MH_PULL_SQUARE_HPP) $(CONFIG_MH_PULL_CUBIC_HPP) $(CONFIG_MH_VSHD_SQUARE_HPP) $(CONFIG_SA_TM_SQUARE_HPP)
