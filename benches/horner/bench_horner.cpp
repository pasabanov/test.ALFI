#include <benchmark/benchmark.h>

#include "../bench_utils.h"

static const std::vector<std::pair<std::vector<double>,std::vector<double>>> polynomials = {
	{{2, -1, 2, -3, 1, -4}, {2, 19, 74, 145, 141, 50}},
	{{1, 1, -1}, {1, 3, 1}},
	{{1, -2, -1, -4, -3, 1, 1, -1}, {1, -16, 107, -394, 869, -1149, 841, -263}},
};

static const std::vector<double> X = {-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
static const std::vector<double> Y = {-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

static std::vector<double> shift_with_binomials(const std::vector<double>& P, const std::vector<double>& X, const std::vector<double>& Y) {
	const auto n = P.size();

	const static auto binomials = [](size_t m) {
		std::vector<std::vector<double>> C(m + 1);
		for (size_t i = 0; i <= m; ++i) {
			C[i].resize(i + 1);
			C[i][0] = C[i][i] = 1;
			for (size_t j = 1; j <= i / 2; ++j) {
				C[i][j] = C[i][i-j] = C[i-1][j-1] + C[i-1][j];
			}
		}
		return C;
	};

	const auto C = binomials(n - 1);

	std::vector<double> coeffs(n * X.size());

	for (size_t i = 0; i < X.size() - 1; ++i) {
		coeffs[i*n] = P[0];
		for (size_t k = 1; k < n - 1; ++k) {
			double r = 0;
			for (size_t j = 0; j < k; ++j) {
				r += ((((k - j) & 1) == 0) ? 1 : -1) * coeffs[i*n+j] * C[n-1-j][k-j];
				r *= X[i];
			}
			coeffs[i*n+k] = P[k] - r;
		}
		coeffs[i*n+n-1] = Y[i];
	}
	return coeffs;
}

static std::vector<double> shift_with_horner(const std::vector<double>& P, const std::vector<double>& X, const std::vector<double>& Y) {
	const auto n = P.size();

	std::vector<double> coeffs(n * X.size());

	for (size_t i = 0; i < X.size() - 1; ++i) {
		for (size_t j = 0; j < n - 1; ++j) {
			coeffs[i*n+j] = P[j];
		}
		for (size_t k = 1; k < n - 1; ++k) {
			for (size_t j = 1; j < n - 1 - k; ++j) {
				coeffs[i*n+j] += coeffs[i*n+j-1] * X[i];
			}
		}
		coeffs[i*n+n-1] = Y[i];
	}
	return coeffs;
}

template <auto shift_function>
static void BM_shift(benchmark::State& state) {
	const auto& [original, expected] = polynomials[state.range(0)];
	std::vector<double> result;
	for (auto _ : state) {
		result = shift_function(original, X, Y);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	std::vector<double> expected_coeffs(expected.size() * X.size());
	for (size_t i = 0; i < X.size(); ++i) {
		std::ranges::copy(expected, expected_coeffs.begin() + static_cast<int64_t>(i * expected.size()));
	}
	add_error_metrics(state, result, expected_coeffs);
}

static void BM_shift_with_binomials(benchmark::State& state) {
	BM_shift<shift_with_binomials>(state);
}

static void BM_shift_with_horner(benchmark::State& state) {
	BM_shift<shift_with_horner>(state);
}

BENCHMARK(BM_shift_with_binomials)
	->ArgsProduct({benchmark::CreateDenseRange(0, static_cast<int64_t>(polynomials.size()) - 1, 1)})
	->ArgName("polynomial");

BENCHMARK(BM_shift_with_horner)
	->ArgsProduct({benchmark::CreateDenseRange(0, static_cast<int64_t>(polynomials.size()) - 1, 1)})
	->ArgName("polynomial");

BENCHMARK_MAIN();