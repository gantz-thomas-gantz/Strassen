## Project Structure

- `src/`: Contains source code files.
- `include/`: Contains header files.
- `CMakeLists.txt`: Defines the build configuration.
- `build/`: Created during building of CMake project and eventually holds the compiled binary.
- `data/`: Data containing test results.
- `README.md`: This file.

## Prerequisites

Ensure you have the following installed on your system:

- GCC or Clang (for compiling C code)
- CMake (version 3.10 or higher recommended)
- CBLAS and LAPACK. Ubunto: sudo apt install libopenblas-dev liblapack-dev


## Compilation

Follow these steps to compile and run the project:

1. Go to the build directory of the project.
   ```bash
   cd build
   ```

2. Run the compilation of the project:
   ```bash
   make
   ```

3. Run the executable to run the tests:
   ```bash
   ./main
   ```

## Notes

- If you want to enable optimizations or see warnings, the project already configures them by default:
  - Optimization: Level 3 (`-O3`)
  - All warnings: (`-Wall -Wextra`)




