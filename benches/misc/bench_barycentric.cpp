#include <benchmark/benchmark.h>

#include <ALFI/misc.h>
#include <ALFI/util/arrays.h>
#include <ALFI/util/stat.h>

#include "../bench_utils.h"
#include "../bench_data.h"

static void BM_barycentric(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const auto& nn = state.range(3);
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	const std::vector<double> xx = of_type<double>(dist_type, nn, interval.first, interval.second);
	std::vector<double> result;
	for (auto _ : state) {
		result = alfi::misc::barycentric(X, Y, xx, dist_type);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	const std::vector<double> yy = apply_func(xx, lambda);
	const std::vector<double> error = alfi::util::arrays::sub(result, yy);
	state.counters["1.ME"] = alfi::util::stat::mean(error);
	state.counters["2.MAE"] = alfi::util::stat::mean_abs(error);
	state.counters["3.RMSE"] = alfi::util::stat::rms(error);
	state.counters["4.Variance"] = alfi::util::stat::var(error);
}
BENCHMARK(BM_barycentric)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		benchmark::CreateRange(8, 32, 2),
		{nn}
	})
->ArgNames({"func", "dist", "n", "nn"});

BENCHMARK(BM_barycentric)
	->Name("BM_barycentric_chebyshev")
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		{1},
		benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3),
		{nn}
	})
->ArgNames({"func", "dist", "n", "nn"});

BENCHMARK(BM_barycentric)
	->Name("BM_barycentric_chebyshev_2")
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		{2},
		benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3),
		{nn}
	})
	->ArgNames({"func", "dist", "n", "nn"});

BENCHMARK_MAIN();