# Compiler and tools
CC = /usr/bin/cc
RM = /bin/rm
CFLAGS = -O2
DEBUG_FLAGS = -O0 -g  # Debug flags with no optimization and debug symbols

# Define directories for build and binary files
BUILD_DIR = build
BIN_DIR = bin

# Object files in the build directory
LIBRARY = $(BUILD_DIR)/matrix.o $(BUILD_DIR)/L2_distance.o $(BUILD_DIR)/matrixadv.o $(BUILD_DIR)/qr.o $(BUILD_DIR)/eigen.o $(BUILD_DIR)/qsort.o $(BUILD_DIR)/svd.o

# Test executables in the bin directory
TEST_APS = $(BIN_DIR)/qr_test $(BIN_DIR)/invtest $(BIN_DIR)/eigen_test $(BIN_DIR)/quicksort_test $(BIN_DIR)/svd_test

# Default target: all
all: $(BUILD_DIR) $(BIN_DIR) $(TEST_APS)

# Debug target: builds with debug flags
debug: CFLAGS += $(DEBUG_FLAGS)
debug: clean all

# Create build and bin directories if they don't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Compile .c files to .o files in the build directory
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile test executables and place them in the bin directory
$(BIN_DIR)/%: %.c $(LIBRARY) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(LIBRARY) $< -o $@ -lm -lgsl -lgslcblas

# Build a specific target if TARGET is provided
ifneq ($(TARGET),)
TARGET_EXEC = $(BIN_DIR)/$(TARGET)

$(TARGET_EXEC): $(TARGET).c $(LIBRARY) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(LIBRARY) $< -o $@ -lm

.PHONY: target
target: $(TARGET_EXEC)
endif

# Clean up build and bin directories
clean:
	-$(RM) -r $(BUILD_DIR)
	-$(RM) -r $(BIN_DIR)
