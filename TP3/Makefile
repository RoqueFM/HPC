# Directorio de los archivos de cabecera (headers)
IDIR = /usr/include/qt5

# Directorio de salida para los objetos
ODIR = .

# Bibliotecas necesarias para enlazar
LIBS = -L/usr/lib/x86_64-linux-gnu #-lQt5Core -lm

# Compilador
CXX = g++

# Flags del compilador
CXXFLAGS = -I$(IDIR) -Wall -fopenmp -O3

# Nombre del ejecutable final
TARGET = main

# Archivos fuente
SRCS = main.cpp

# Archivos objeto
OBJS = $(SRCS:.cpp=.o)

# Regla por defecto: compilar y ejecutar
all: $(TARGET)

# Regla para enlazar los archivos objeto en el ejecutable final
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Regla para compilar archivos .cpp en archivos .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Regla para ejecutar el programa
run: $(TARGET)
	./$(TARGET)

# Regla para limpiar archivos generados
clean:
	rm -f $(TARGET) $(OBJS)
