CC = g++
#CC = clang++

CXXFLAGS = -g
LIBS = -lgmp -lgmpxx

SRCS = $(wildcard src/*.cpp)

HEADERS = $(wildcard include/*.h)

TEST_MAIN = tests/main.cpp

INPUT_FILES = $(wildcard tests/input/*)

OUTPUT = a.out

OBJS = $(SRCS:.cpp=.o)

all: $(OUTPUT)

$(OUTPUT): $(OBJS) $(TEST_MAIN)
	$(CC) $(CXXFLAGS) $(OBJS) $(TEST_MAIN) -o $(OUTPUT) $(LIBS)

src/%.o: src/%.cpp $(HEADERS)
	$(CC) $(CXXFLAGS) -c $< -o $@

run: $(OUTPUT)
	@for input_file in $(INPUT_FILES); do \
		echo "#################"; \
		echo $$input_file; \
		echo "#################"; \
		echo "expected output:"; \
		output_file=tests/output/$$(basename $$input_file); \
		cat $$output_file; \
		echo "computed output:"; \
		./$(OUTPUT) $$input_file tmp.output; \
		cat tmp.output; \
		cmp --silent $$output_file tmp.output  && echo '### PASSED ###' || echo '### FAILED! ###'; \
		rm -f tmp.output; \
		echo ""; \
	done

clean:
	rm -f src/*.o $(OUTPUT)

.PHONY: clean all

