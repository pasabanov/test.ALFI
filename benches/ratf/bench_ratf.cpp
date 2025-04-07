#include <benchmark/benchmark.h>

#include <ALFI/ratf.h>

#include "../bench_utils.h"

static const std::vector<std::pair<double,double>> intervals = {
	{-1, 1},
	{5, 10},
	{0, 2},
};

static void BM_val_scalar(benchmark::State& state) {
	const auto& n = state.range(0);
	const auto& interval = intervals[state.range(1)];
	alfi::ratf::RationalFunction<> rf;
	rf.first.resize(n);
	rf.second.resize(n);
	for (auto& c : rf.first) {
		c = random_value<double>(-1, 1);
	}
	for (auto& c : rf.second) {
		c = random_value<double>(-1, 1);
	}
	const auto x = random_value<double>(interval.first, interval.second);
	double result;
	for (auto _ : state) {
		result = alfi::ratf::val(rf, x);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_val_scalar)
	->ArgsProduct({
		benchmark::CreateRange(1 << 3, 1 << 12, 1 << 3),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(intervals.size()) - 1, 1),
	})
	->ArgNames({"n", "interval"});

static void BM_val_vector(benchmark::State& state) {
	const auto& n = state.range(0);
	const auto& nn = state.range(1);
	const auto& interval = intervals[state.range(2)];
	alfi::ratf::RationalFunction<> rf;
	rf.first.resize(n);
	rf.second.resize(n);
	for (auto& c : rf.first) {
		c = random_value<double>(-1, 1);
	}
	for (auto& c : rf.second) {
		c = random_value<double>(-1, 1);
	}
	std::vector<double> xx(nn);
	for (auto& x : xx) {
		x = random_value<double>(interval.first, interval.second);
	}
	std::vector<double> result;
	for (auto _ : state) {
		result = alfi::ratf::val(rf, xx);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_val_vector)
	->ArgsProduct({
		benchmark::CreateRange(1 << 3, 1 << 12, 1 << 3),
		benchmark::CreateRange(1 << 3, 1 << 12, 1 << 3),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(intervals.size()) - 1, 1),
	})
	->ArgNames({"n", "nn", "interval"});

static void BM_pade(benchmark::State& state) {
	const auto& size = state.range(0);
	const auto degree = size - 1;
	const auto n = degree/2, m = n;
	std::vector<double> p(size);
	for (size_t i = 0; i < p.size(); ++i) {
		p[i] = random_value<double>(-2, 2) * static_cast<double>(i + 1);
	}
	alfi::ratf::RationalFunction<> result;
	for (auto _ : state) {
		result = alfi::ratf::pade(p, n, m);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	state.SetBytesProcessed(state.iterations() * size * static_cast<int64_t>(sizeof(double)));
	state.SetComplexityN(size);
}
BENCHMARK(BM_pade)
	->ArgsProduct({benchmark::CreateRange(1, 1 << 16, 2)})
	->ArgName("n")
	->Complexity();

BENCHMARK_MAIN();