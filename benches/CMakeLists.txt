set(CMAKE_CXX_STANDARD 23)

set(ORIGINAL_CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
set(CMAKE_BUILD_TYPE Release)
set(BENCHMARK_ENABLE_INSTALL OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_LTO ON)
add_subdirectory(deps/benchmark)
set(CMAKE_BUILD_TYPE "${ORIGINAL_CMAKE_BUILD_TYPE}")

add_executable(benchmarks
	poly/bench_poly.cpp)
target_link_libraries(benchmarks PRIVATE ALFI benchmark::benchmark)

add_custom_target(bench
	COMMAND benchmarks --benchmark_color=yes --benchmark_counters_tabular=yes
	DEPENDS benchmarks
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Running benchmarks..."
	USES_TERMINAL)