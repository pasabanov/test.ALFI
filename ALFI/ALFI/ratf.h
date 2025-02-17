#pragma once

#include <iostream>
#include <cmath>

#include "config.h"
#include "util/numeric.h"

namespace alfi::ratf {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	using RationalFunction = std::pair<Container<Number>, Container<Number>>;

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val_mul(const RationalFunction<Number, Container>& rf, Number x) {
		Number n = 0;
		for (const Number& c : rf.first) {
			n = x * n + c;
		}
		Number d = 0;
		for (const Number& c : rf.second) {
			d = x * d + c;
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
	RationalFunction<Number,Container> pade(const Container<Number> P, SizeT n, SizeT m) {
		// const auto get_p = [&](SizeT i) {
		// 	if (i >= P.size()) {
		// 		return 0;
		// 	}
		// 	return P[i];
		// };
		// const auto d = P.size();
		// Container<Number> A(n);
		// Container<Number> B(m);
		// // evaluating B
		//
		// // evaluating A
		// for (SizeT i = 0; i < n; ++i) {
		// 	A[i] = get_p(i);
		// 	for (SizeT j = 0; j <= i; ++j) {
		// 		A[i] += B[j] * get_p(i-1-j);
		// 	}
		// }
		// return {A, B};
	}
}