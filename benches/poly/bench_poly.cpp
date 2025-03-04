#include <cmath>
#include <limits>
#include <vector>

#include <benchmark/benchmark.h>

#include "ALFI/poly.h"

static void BM_val_scalar(benchmark::State& state) {
	const std::vector<double> coeffs{1, 2, 3, 4, 5};
	const double x = 1.234;
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::val(coeffs, x));
	}
	state.counters["abc"] = 1;
	state.counters["abcd"] = 1;
	state.counters["qwe"] = 2;
}
BENCHMARK(BM_val_scalar);

static void BM_val_vector(benchmark::State& state) {
	const std::vector<double> coeffs{1, 2, 3, 4, 5};
	std::vector<double> xx(1000);
	for (std::size_t i = 0; i < xx.size(); ++i) {
		xx[i] = 1.0 + i * 0.01;
	}
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::val(coeffs, xx));
	}
	state.counters["abc"] = 1;
	state.counters["qwe"] = 2;
}
BENCHMARK(BM_val_vector);

static void BM_lagrange(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::lagrange(X, Y));
	}
}
BENCHMARK(BM_lagrange);

static void BM_lagrange_vals(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	std::vector<double> xx(1000);
	for (std::size_t i = 0; i < xx.size(); ++i) {
		xx[i] = 0.0 + i * 0.009;
	}
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::lagrange_vals(X, Y, xx));
	}
}
BENCHMARK(BM_lagrange_vals);

static void BM_imp_lagrange(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::imp_lagrange(X, Y));
	}
}
BENCHMARK(BM_imp_lagrange);

static void BM_imp_lagrange_vals(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	std::vector<double> xx(1000);
	for (std::size_t i = 0; i < xx.size(); ++i) {
		xx[i] = 0.0 + i * 0.009;
	}
	double epsilon = std::numeric_limits<double>::epsilon();
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::imp_lagrange_vals(X, Y, xx, epsilon));
	}
}
BENCHMARK(BM_imp_lagrange_vals);

static void BM_newton(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::newton(X, Y));
	}
}
BENCHMARK(BM_newton);

static void BM_newton_vals(benchmark::State& state) {
	const std::vector<double> X{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	const std::vector<double> Y{0, 1, 4, 9, 16, 25, 36, 49, 64, 81};
	std::vector<double> xx(1000);
	for (std::size_t i = 0; i < xx.size(); ++i) {
		xx[i] = 0.0 + i * 0.009;
	}
	for (auto _ : state) {
		benchmark::DoNotOptimize(alfi::poly::newton_vals(X, Y, xx));
	}
}
BENCHMARK(BM_newton_vals);

BENCHMARK_MAIN();