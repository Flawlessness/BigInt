add_executable(BigIntTest main.cpp BigIntTest.cpp)

target_link_libraries(BigIntTest PRIVATE BigIntLib gtest)
