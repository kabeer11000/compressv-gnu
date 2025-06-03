# Compiler
CXX = riscv64-linux-gnu-g++


# Path to the directory that contains the 'evpp' header root
# This is where 'evpp/tcp_server.h' is found relative to.
# Since your main directory is 'pod', and 3rdparty is inside 'pod',
# this path correctly points to 'pod/3rdparty/evpp'.

# Compiler flags
# -I$(EVPP_INCLUDE_ROOT) tells g++ to look in '3rdparty/evpp' for headers.
# So, '#include <evpp/tcp_server.h>' will correctly resolve to '3rdparty/evpp/evpp/tcp_server.h'.
CXXFLAGS = -S -O0 -march=rv64gcv -mabi=lp64d -static -pthread -lrt -I thirdparty/crow/include

# Linker flags
# These link your executable with the evpp library and its dependencies.
LDFLAGS = 

# Target executable
TARGET = build/compressv

# For deleting the target and object files
TARGET_DEL = $(TARGET) $(OBJS)

# Source files
# Assuming server.cpp is directly inside the 'pod' main directory.
SRCS = $(wildcard compressv/*.c)

# Object files
OBJS = $(SRCS:.c=.o)

# Default rule to build and run the executable
all: $(TARGET) run

# Rule to link object files into the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Rule to compile .cpp files into .o files
%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to run the executable
run: $(TARGET)
	./$(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(TARGET_DEL)