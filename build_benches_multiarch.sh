cmake \
	-S . \
	-B build/x86_64-debug \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_CXX_COMPILER=x86_64-linux-gnu-g++ \
	-DCMAKE_OBJDUMP=x86_64-linux-gnu-objdump \
	-DCMAKE_EXE_LINKER_FLAGS="-static" \
	-DALFI_ENABLE_EXAMPLES=OFF \
	-DALFI_ENABLE_TESTS=OFF \
	-DALFI_ENABLE_BENCHES=ON

cmake --build build/x86_64-debug -j

cmake \
	-S . \
	-B build/riscv64-debug \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_CXX_COMPILER=riscv64-linux-gnu-g++ \
	-DCMAKE_OBJDUMP=riscv64-linux-gnu-objdump \
	-DCMAKE_EXE_LINKER_FLAGS="-static" \
	-DALFI_ENABLE_EXAMPLES=OFF \
	-DALFI_ENABLE_TESTS=OFF \
	-DALFI_ENABLE_BENCHES=ON

cmake --build build/riscv64-debug -j