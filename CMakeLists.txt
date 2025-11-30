cmake_minimum_required(VERSION 3.22)

project(ALFI
	VERSION 0.0.0.0
	DESCRIPTION "Advanced Library for Function Interpolation"
	HOMEPAGE_URL "https://alfi-lib.github.io/ALFI/"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-g3 -O3 -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-g3 -Oz -flto=auto -DNDEBUG")
set(CMAKE_CXX_FLAGS_FAST "-g3 -Ofast -fno-finite-math-only -flto=auto -DNDEBUG")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ALFI_ENABLE_EXAMPLES "Enable examples subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})
option(ALFI_ENABLE_TESTS "Enable tests subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})
option(ALFI_ENABLE_BENCHES "Enable benches subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})

add_subdirectory(ALFI)

if(ALFI_ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()

if(ALFI_ENABLE_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()

if(ALFI_ENABLE_BENCHES)
	add_subdirectory(benches)
endif()