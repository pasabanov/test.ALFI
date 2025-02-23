#pragma once

#include <iostream>
#include <cmath>

#include "config.h"
#include "util/linalg.h"
#include "util/numeric.h"

namespace alfi::ratf {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	using RationalFunction = std::pair<Container<Number>, Container<Number>>;

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val_mul(const RationalFunction<Number, Container>& rf, Number x) {
		Number n = 0;
		for (const auto& c : rf.first) {
			n = n * x + c;
		}
		Number d = 0;
		for (const auto& c : rf.second) {
			d = d * x + c;
		}
		return n / d;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val_mul(const RationalFunction<Number, Container>& rf, const Container<Number>& X) {
		Container<Number> result(X.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < X.size(); ++i) {
			result[i] = val_mul(rf, X[i]);
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val_div(const RationalFunction<Number, Container>& rf, Number x) {
		const auto& numerator = rf.first;
		const auto& denominator = rf.second;
		Number n = 0;
		for (auto i = numerator.rbegin(); i != numerator.rend(); ++i) {
			n = n / x + *i;
		}
		Number d = 0;
		for (auto i = denominator.rbegin(); i != denominator.rend(); ++i) {
			d = d / x + *i;
		}
		const auto numerator_degree = numerator.empty() ? 0 : numerator.size() - 1;
		const auto denominator_degree = denominator.empty() ? 0 : denominator.size() - 1;
		if (numerator_degree >= denominator_degree) {
			return n / d * util::numeric::binpow(x, numerator_degree - denominator_degree);
		} else {
			return n / d / util::numeric::binpow(x, denominator_degree - numerator_degree);
		}
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val_div(const RationalFunction<Number, Container>& rf, const Container<Number>& X) {
		Container<Number> result(X.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < X.size(); ++i) {
			result[i] = val_div(rf, X[i]);
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val(const RationalFunction<Number, Container>& rf, Number x) {
		if (std::abs(x) <= 1) {
			return val_mul(rf, x);
		} else {
			return val_div(rf, x);
		}
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val(const RationalFunction<Number, Container>& rf, const Container<Number>& X) {
		Container<Number> result(X.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < X.size(); ++i) {
			result[i] = val(rf, X[i]);
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	RationalFunction<Number,Container> pade(const Container<Number>& P, SizeT n, SizeT m, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		if (m == 0) {
			return {P, {static_cast<Number>(1)}};
		}
		const static auto get_or_0 = [](const Container<Number>& arr, SizeT index) -> Number {
			return 0 <= index && index < arr.size() ? arr[index] : static_cast<Number>(0);
		};
		Container<Number> B(m + 1);
		// I. Evaluating B
		// Matrix
		Container<Container<Number>> M(m);
		// Right-hand side
		Container<Number> R(m);
		// Filling system of equations
		for (SizeT row = 0; row < m; ++row) {
			R[row] = -get_or_0(P, n + 1 + row);
			for (SizeT col = 0; col < m; ++col) {
				M[row][col] = get_or_0(P, n + row - col);
			}
		}
		// Removing zeroed rows
		for (SizeT row = n + m - 1; row > n; ++row) {
			bool all_zero = true;
			if (std::abs(R[row]) < epsilon) {
				for (SizeT col = 0; col < m; ++col) {
					if (std::abs(M[row][col]) >= epsilon) {
						all_zero = false;
						break;
					}
				}
			}
			if (!all_zero) {
				break;
			}
			B[row+1-n] = 0;
			M.pop_back();
			for (auto& r : M) {
				r.pop_back();
			}
			R.pop_back();
		}
		// Solving system of equations
		auto X = util::linalg::lup_solve(M, R, epsilon);
		if (X.empty()) {
			std::cout << "Failed to construct Pade approximation" << std::endl;
			return {{}, {}};
		}
		std::move(X.begin(), X.end(), B.begin() + 1);
		Container<Number> A(n + 1);
		// II. Evaluating A
		for (SizeT i = 0; i < A.size(); ++i) {
			A[i] = get_or_0(P, i);
			for (SizeT j = 1; j <= i; ++j) {
				A[i] += get_or_0(B, j - 1) * get_or_0(P, i - j);
			}
		}
		return {A, B};
	}
}