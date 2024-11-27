cmake_minimum_required(VERSION 3.12)
project(EX1)
add_executable(main src/main.c src/IO.c src/naive_matmat.c src/strassen_matmat.c src/test.c)
target_include_directories(main PUBLIC include)
target_link_libraries(main -lm)
target_link_libraries(main -lcblas)

