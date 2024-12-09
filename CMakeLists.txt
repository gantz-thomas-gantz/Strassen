cmake_minimum_required(VERSION 3.12)
project(EX1)
add_executable(main src/main.c src/IO.c src/block_utilities.c src/naive_matmat.c 
	src/strassen_matmat.c src/strassen_inv.c src/test.c)
target_include_directories(main PUBLIC include)
target_link_libraries(main PRIVATE m)
target_link_libraries(main PRIVATE lapacke cblas m)
