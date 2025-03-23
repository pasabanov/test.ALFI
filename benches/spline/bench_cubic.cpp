#include <benchmark/benchmark.h>

#include <ALFI/spline/cubic.h>
#include <ALFI/util/arrays.h>
#include <ALFI/util/stat.h>

#include "../bench_utils.h"
#include "../bench_data.h"

static void BM_cubic_spline(benchmark::State& state) {
	for (auto _ : state) {

	}
}

BENCHMARK(BM_cubic_spline)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		benchmark::CreateRange(8, 32, 2),
		{nn}
	})
	->ArgNames({"func", "dist", "n", "nn"});

BENCHMARK_MAIN();