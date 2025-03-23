#include <cmath>
#include <vector>

#include <benchmark/benchmark.h>

#include <ALFI/poly.h>
#include <ALFI/util/arrays.h>
#include <ALFI/util/stat.h>

#include "../bench_utils.h"
#include "../bench_data.h"

static void BM_val_scalar(benchmark::State& state) {
	std::vector<double> coeffs(state.range(0));
	for (auto& coeff : coeffs) {
		coeff = random_value<double>(-1, 1);
	}
	const auto x = random_value<double>(-1, 1);
	double result;
	for (auto _ : state) {
		result = alfi::poly::val(coeffs, x);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
}

static void BM_val_vector(benchmark::State& state) {
	std::vector<double> coeffs(state.range(0));
	for (auto& coeff : coeffs) {
		coeff = random_value<double>(-1, 1);
	}
	std::vector<double> xx(state.range(1));
	for (auto& x : xx) {
		x = random_value<double>(-1, 1);
	}
	std::vector<double> result;
	for (auto _ : state) {
		result = alfi::poly::val(coeffs, xx);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
}

template <auto interp>
static void BM_interp_poly(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const auto& nn = state.range(3);
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	std::vector<double> coeffs;
	for (auto _ : state) {
		coeffs = interp(X, Y);
		benchmark::DoNotOptimize(coeffs);
		benchmark::ClobberMemory();
	}
	const std::vector<double> xx = of_type<double>(dist_type, nn, interval.first, interval.second);
	const std::vector<double> yy = apply_func(xx, lambda);
	const std::vector<double> result = alfi::poly::val(coeffs, xx);
	const std::vector<double> error = alfi::util::arrays::sub(result, yy);
	state.counters["1.ME"] = alfi::util::stat::mean(error);
	state.counters["2.MAE"] = alfi::util::stat::mean_abs(error);
	state.counters["3.RMSE"] = alfi::util::stat::rms(error);
	state.counters["4.Variance"] = alfi::util::stat::var(error);
}

template <auto interp>
static void BM_interp_poly_vals(benchmark::State& state) {
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
		result = interp(X, Y, xx);
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

// need special function for `imp_lagrange_vals` because of the fourth parameter
static void BM_imp_lagrange_vals(benchmark::State& state) {
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
		result = alfi::poly::imp_lagrange_vals(X, Y, xx);
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

struct BenchmarkRegistrar {
	BenchmarkRegistrar() {
		RegisterBenchmark("BM_val_scalar", BM_val_scalar)->Range(1 << 3, 1 << 12);
		RegisterBenchmark("BM_val_vector", BM_val_vector)->Ranges({{1 << 6, 1 << 12}, {1 << 6, 1 << 12}});

		const std::vector<std::pair<std::string, std::function<void(benchmark::State&)>>> regsPoly = {
			{"BM_lagrange", [](benchmark::State& state) { BM_interp_poly<alfi::poly::lagrange<>>(state); }},
			{"BM_imp_lagrange", [](benchmark::State& state) { BM_interp_poly<alfi::poly::imp_lagrange<>>(state); }},
			{"BM_newton", [](benchmark::State& state) { BM_interp_poly<alfi::poly::newton<>>(state); }},
			{"BM_lagrange_vals", [](benchmark::State& state) { BM_interp_poly_vals<alfi::poly::lagrange_vals<>>(state); }},
			{"BM_imp_lagrange_vals", BM_imp_lagrange_vals},
			{"BM_newton_vals", [](benchmark::State& state) { BM_interp_poly_vals<alfi::poly::newton_vals<>>(state); }},
		};
		for (const auto& [name, function] : regsPoly) {
			RegisterBenchmark(name, function)
				->ArgsProduct({
					benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
					benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
					benchmark::CreateRange(8, 32, 2),
					{nn}
				})
				->ArgNames({"func", "dist", "n", "nn"});
		}
	}
};

[[maybe_unused]] const BenchmarkRegistrar registrar;