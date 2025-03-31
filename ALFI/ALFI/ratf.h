#pragma once

#include <cmath>

#include "config.h"
#include "util/numeric.h"
#include "util/poly.h"

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

	/**
		@brief Computes the Padé approximation of a given power series.

		This function constructs a Padé approximant \f(R(x) = P_n(x) / Q_m(x)\f) for a given power series
		represented by its coefficients \f(P\f). The approximant is determined by finding polynomials
		\f(P_n\f) and \f(Q_m\f) of degrees at most \f(n\f) and \f(m\f), respectively, such that:

		\f[
			Q_m(x) \cdot \sum_{k=0}^{n+m} P_k x^k - P_n(x) = \mathcal{O}(x^{n+m+1})
		\f]

		To achieve this, the algorithm solves a linear system using the extended Euclidean algorithm for polynomials.
		The resulting polynomials \f(P_n\f) and \f(Q_m\f) are guaranteed to be normalized (i.e., leading coefficients are nonzero).

		If no valid Padé approximation can be found (e.g., due to a degenerate system), the function returns empty polynomials.

		@param P The power series coefficients in descending order (from highest to lowest degree).
		@param n The desired degree of the numerator polynomial \f(P_n\f).
		@param m The desired degree of the denominator polynomial \f(Q_m\f).
		@param epsilon A threshold for numerical stability when determining polynomial degree (default: machine epsilon).
		@return A rational function \f(P_n(x) / Q_m(x)\f), or an empty result \f({}, {}\f) if the approximation is not possible.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	RationalFunction<Number,Container> pade(const Container<Number>& P, SizeT n, SizeT m, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		if (P.size() < n + m + 1) {
			Container<Number> P1 = P;
			P1.insert(P1.begin(), n + m + 1 - P.size(), 0);
			return pade(P1, n, m, epsilon);
		}

		Container<Number> Xmn1(n + m + 2, 0);
		Xmn1[0] = 1;

		auto [r, _, t] = util::poly::extended_euclid(Xmn1, P, epsilon, n);

		if (r.size() > n + 1 || t.size() > m + 1) {
			return {{}, {}};
		}

		return RationalFunction<Number,Container>(std::move(r), std::move(t));
	}
}