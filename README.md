AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen

## Project Structure

- `src/`: Contains source code files.
- `include/`: Contains header files.
- `CMakeLists.txt`: Defines the build configuration.
- `build/`: Created during building of CMake project and eventually holds the compiled binary.
- `README.md`: This file.

## Prerequisites

Ensure you have the following installed on your system:

- GCC or Clang (for compiling C code)
- CMake (version 3.10 or higher recommended)
- CBLAS and LAPACK. Ubuntu: sudo apt install libopenblas-dev liblapack-dev

## Compilation

Follow these steps to compile and run the project:

1. Run the compilation of the project:
   ```bash
   cmake -B build
   cmake --build build
   ```

2. Run the executable to run the tests:
   ```bash
   ./main <test size (default 5)>
   ```

3. See results in console or optionally in build/matinv.txt and
   build/matmat.txt.

## Notes

- If you want to enable optimizations or see warnings, the project already configures them by default:
  - Optimization: Level 3 (`-O3`)
  - All warnings: (`-Wall -Wextra`)




