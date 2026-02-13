CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

TARGET = main
SRCS = main.cpp int512.cpp float512.cpp
OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
