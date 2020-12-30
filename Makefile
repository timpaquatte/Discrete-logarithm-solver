SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

EXE := $(BIN_DIR)/solve_dlog
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)


LD_FLAGS=-lgmp -lm
CC_FLAGS=-O3
INC=-Iinclude

all: $(EXE)

$(EXE):	$(OBJ) | $(BIN_DIR)
	$(CC) $^ -o $@ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(INC) -c $^ $(CC_FLAGS) -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	rm -rv $(BIN_DIR) $(OB_DIR)