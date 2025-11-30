set(CMAKE_CXX_STANDARD 23)

set(ORIGINAL_CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
set(CMAKE_BUILD_TYPE Release)
set(BENCHMARK_ENABLE_INSTALL OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_LTO ON)
set(HAVE_STD_REGEX TRUE) # Hack for RISC-V
add_subdirectory(deps/benchmark)
set(CMAKE_BUILD_TYPE "${ORIGINAL_CMAKE_BUILD_TYPE}")

macro(add_bench_executable bench_name)
	add_executable(${bench_name} ${ARGN} bench_utils.h bench_data.h)
	target_compile_options(${bench_name} PRIVATE -fverbose-asm -save-temps=obj)
	target_link_libraries(${bench_name} PRIVATE ALFI benchmark::benchmark)
	add_custom_command(
		TARGET ${bench_name}
		POST_BUILD
		COMMAND ${CMAKE_OBJDUMP} -d -S -C $<TARGET_FILE:${bench_name}> > $<TARGET_FILE_DIR:${bench_name}>/${bench_name}.objdump.s
		VERBATIM
	)
endmacro()

add_bench_executable(bench_poly poly/bench_poly.cpp)
add_bench_executable(bench_ratf ratf/bench_ratf.cpp)
add_bench_executable(bench_barycentric misc/bench_barycentric.cpp)
add_bench_executable(bench_polyeqv spline/bench_polyeqv.cpp)
add_bench_executable(bench_step spline/bench_step.cpp)
add_bench_executable(bench_linear spline/bench_linear.cpp)
add_bench_executable(bench_quadratic spline/bench_quadratic.cpp)
add_bench_executable(bench_cubic spline/bench_cubic.cpp)