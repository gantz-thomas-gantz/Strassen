cmake_minimum_required(VERSION 3.12)
project(EX1)
add_executable(main src/main.c src/IO.c src/block_utilities.c src/naive_matmat.c 
	src/strassen_matmat.c src/strassen_inv.c src/naive_lu.c src/test.c)
target_include_directories(main PUBLIC include)
target_link_libraries(main PRIVATE m)
target_link_libraries(main PRIVATE lapacke cblas m)

# Set optimization level to 3
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# Enable all warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

