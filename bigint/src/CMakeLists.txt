add_executable(BigInt main.cpp)

add_library(BigIntLib BigInt.cpp BigInt.h)

target_link_libraries(BigInt PRIVATE BigIntLib)