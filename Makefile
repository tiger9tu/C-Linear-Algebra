CC = /usr/bin/cc
RM = /bin/rm
CFLAGS = -O2

# Define directories for build and binary files
BUILD_DIR = build
BIN_DIR = bin

# Object files in the build directory
LIBRARY = $(BUILD_DIR)/matrix.o $(BUILD_DIR)/L2_distance.o $(BUILD_DIR)/matrixadv.o $(BUILD_DIR)/qr.o $(BUILD_DIR)/eigen.o $(BUILD_DIR)/qsort.o $(BUILD_DIR)/svd.o

# Test executables in the bin directory
TEST_APS = $(BIN_DIR)/qr_decomposition_test $(BIN_DIR)/invtest $(BIN_DIR)/eigen_test $(BIN_DIR)/quicksort_test $(BIN_DIR)/svd_test

# Target: all
all: $(BUILD_DIR) $(BIN_DIR) $(TEST_APS)

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
	$(CC) $(CFLAGS) $(LIBRARY) $< -o $@ -lm

# Clean up build and bin directories
clean:
	-$(RM) -r $(BUILD_DIR)
	-$(RM) -r $(BIN_DIR)
