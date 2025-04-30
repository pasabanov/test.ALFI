#include <benchmark/benchmark.h>

#include "../bench_utils.h"

static const std::vector<std::tuple<std::vector<double>,double,std::vector<double>>> polynomials = {
	// x^2 + x - 1 by (x-1)
	{{1, 1, -1}, 1, {1, 3, 1}},
	// 2x^5 -x^4 + 2x^3 - 3x^2 + x - 4 by (x-2)
	{{2, -1, 2, -3, 1, -4}, 2, {2, 19, 74, 145, 141, 50}},
	// x^7 - 2x^6 - x^5 - 4x^4 - 3x^3 + x^2 + x - 1 by (x+2)
	{{1, -2, -1, -4, -3, 1, 1, -1}, -2, {1, -16, 107, -394, 869, -1149, 841, -263}},
	// x^20 - 1 by (x-1)
	{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1}, 1, {1,-20,190,-1140,4845,-15504,38760,-77520,125970,-167960,184756,-167960,125970,-77520,38760,-15504,4845,-1140,190,-20, 0}},
	// x^30 - 1 by (x-1)
	{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1}, 1, {1, -30, 435, -4060, 27405, -142506, 593775, -2035800, 5852925, -14307150, 30045015, -54627300, 86493225, -119759850, 145422675, -155117520, 145422675, -119759850, 86493225, -54627300, 30045015, -14307150, 5852925, -2035800, 593775, -142506, 27405, -4060, 435, -30, 0}},
	// x^40 - 1 by (x-1)
	{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1}, 1, {1, -40, 780, -9880, 91390, -658008, 3838380, -18643560, 76904685, -273438880, 847660528, -2311801440, 5586853480, -12033222880, 23206929840, -40225345056, 62852101650, -88732378800, 113380261800, -131282408400, 137846528820, -131282408400, 113380261800, -88732378800, 62852101650, -40225345056, 23206929840, -12033222880, 5586853480, -2311801440, 847660528, -273438880, 76904685, -18643560, 3838380, -658008, 91390, -9880, 780, -40, 0}},
};

static const double Y = 5; // some number

static std::vector<std::vector<double>> binomials(size_t m) {
	std::vector<std::vector<double>> C(m + 1);
	for (size_t i = 0; i <= m; ++i) {
		C[i].resize(i + 1);
		C[i][0] = C[i][i] = 1;
		for (size_t j = 1; j <= i / 2; ++j) {
			C[i][j] = C[i][i-j] = C[i-1][j-1] + C[i-1][j];
		}
	}
	return C;
}

static std::vector<double> shift_with_binomials(const std::vector<double>& P, double x0) {
	const auto n = P.size();
	const auto C = binomials(n - 1);
	std::vector<double> result(n);

	result[0] = P[0];
	for (size_t k = 1; k < n - 1; ++k) {
		double r = 0;
		for (size_t j = 0; j < k; ++j) {
			r += ((((k - j) & 1) == 0) ? 1 : -1) * result[j] * C[n-1-j][k-j];
			r *= x0;
		}
		result[k] = P[k] - r;
	}
	result[n-1] = Y;

	return result;
}

static std::vector<double> shift_with_binomials_C_predefined(const std::vector<double>& P, double x0, const std::vector<std::vector<double>>& C) {
	const auto n = P.size();
	std::vector<double> result(n);

	result[0] = P[0];
	for (size_t k = 1; k < n - 1; ++k) {
		double r = 0;
		for (size_t j = 0; j < k; ++j) {
			r += ((((k - j) & 1) == 0) ? 1 : -1) * result[j] * C[n-1-j][k-j];
			r *= x0;
		}
		result[k] = P[k] - r;
	}
	result[n-1] = Y;

	return result;
}

static std::vector<double> shift_with_horner(const std::vector<double>& P, double x0) {
	const auto n = P.size();
	std::vector<double> result(n);

	for (size_t j = 0; j < n; ++j) {
		result[j] = P[j];
	}
	for (size_t k = 0; k < n - 1; ++k) {
		for (size_t j = 1; j < n - k; ++j) {
			result[j] += result[j-1] * x0;
		}
	}
	result[n-1] = Y;

	return result;
}

static void BM_shift_with_binomials(benchmark::State& state) {
	const auto& [original, pivot, expected] = polynomials[state.range(0)];
	std::vector<double> result;
	for (auto _ : state) {
		result = shift_with_binomials(original, pivot);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	add_error_metrics(state, result, expected);
}

static void BM_shift_with_binomials_C_predefined(benchmark::State& state) {
	const auto& [original, pivot, expected] = polynomials[state.range(0)];
	const auto C = binomials(original.size() - 1);
	std::vector<double> result;
	for (auto _ : state) {
		result = shift_with_binomials_C_predefined(original, pivot, C);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	add_error_metrics(state, result, expected);
}

static void BM_shift_with_horner(benchmark::State& state) {
	const auto& [original, pivot, expected] = polynomials[state.range(0)];
	std::vector<double> result;
	for (auto _ : state) {
		result = shift_with_horner(original, pivot);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	add_error_metrics(state, result, expected);
}

BENCHMARK(BM_shift_with_binomials)
	->ArgsProduct({benchmark::CreateDenseRange(0, static_cast<int64_t>(polynomials.size()) - 1, 1)})
	->ArgName("polynomial");

BENCHMARK(BM_shift_with_binomials_C_predefined)
	->ArgsProduct({benchmark::CreateDenseRange(0, static_cast<int64_t>(polynomials.size()) - 1, 1)})
	->ArgName("polynomial");

BENCHMARK(BM_shift_with_horner)
	->ArgsProduct({benchmark::CreateDenseRange(0, static_cast<int64_t>(polynomials.size()) - 1, 1)})
	->ArgName("polynomial");

BENCHMARK_MAIN();