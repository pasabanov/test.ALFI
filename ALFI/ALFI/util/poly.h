#pragma once

#include <algorithm>
#include <limits>

#include "../config.h"

namespace alfi::util::poly {
	/**
		@brief Normalizes a polynomial by removing leading coefficients close to zero.

		This function removes coefficients from the beginning of the polynomial (i.e. from the highest degree)
		that are smaller than a given threshold \p epsilon. If all coefficients are removed, a single zero is added.

		@param p the polynomial to normalize (a container of coefficients)
		@param epsilon the threshold to treat coefficients as zero (default is machine epsilon)
	 */
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	void normalize(Container<Number>& p, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		if (p.empty()) {
			return p.push_back(0);
		}
		auto p_start = std::find_if(p.begin(), p.end(), [&epsilon](Number v) { return std::abs(v) > epsilon; });
		if (p_start == p.end()) {
			--p_start;
		}
		if (p_start > p.begin()) {
			p.erase(p.begin(), p_start);
		}
	}

	/**
		@brief Multiplies two polynomials.

		Given two polynomials \f(p_1\f) and \f(p_2\f) represented as containers of coefficients,
		this function computes their product using the standard convolution formula:
		\f[
			(p_1 \cdot p_2)[k] = \sum_{i+j=k} p_1[i]\,p_2[j]
		\f]
		where the index corresponds to the degree of the resulting term.

		If either polynomial is empty, the function returns an empty container.

		@param p1 the first polynomial
		@param p2 the second polynomial
		@return the product polynomial \f(p_1 \cdot p_2\f)
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> mul(const Container<Number>& p1, const Container<Number>& p2) {
		if (p1.empty() || p2.empty()) {
			return {};
		}
		Container<Number> result(p1.size() + p2.size() - 1);
		for (SizeT i = 0; i < p1.size(); ++i) {
			for (SizeT j = 0; j < p2.size(); ++j) {
				result[i+j] += p1[i] * p2[j];
			}
		}
		return result;
	}

	/**
		@brief Performs polynomial division.

		This function divides the dividend polynomial \f(A\f) by the divisor polynomial \f(B\f)
		(whose coefficients are in descending order), and returns the quotient \f(Q\f) and the remainder \f(R\f)
		such that:
		\f[
			A(x) = B(X) \cdot Q(x) + R(x)
		\f]
		with \f(deg(R) < deg(B)\f).

		If the divisor is the zero polynomial (i.e. the container is empty or all its coefficients are below \p epsilon)
		or if the dividend has a lower degree than the divisor, the function returns an empty quotient and the original dividend as remainder.

		@param dividend the polynomial to be divided
		@param divisor the polynomial to divide by
		@param epsilon a threshold for detecting zero coefficients (default is machine epsilon)
		@return a pair `{quotient, remainder}`
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::pair<Container<Number>,Container<Number>> div(const Container<Number>& dividend, const Container<Number>& divisor, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		const auto divisor_start
			= std::find_if(divisor.begin(), divisor.end(), [&epsilon](Number v) { return std::abs(v) > epsilon; });

		if (divisor_start == divisor.end() || dividend.size() < divisor.size()) {
			return {{}, dividend};
		}

		const auto divisor_start_idx = divisor_start - divisor.begin();

		const auto n = dividend.size();
		const auto m = divisor.size() - divisor_start_idx;

		Container<Number> quotient(n - m + 1, 0);
		Container<Number> remainder = dividend;

		for (SizeT i = 0; i <= n - m; ++i) {
			const Number factor = remainder[i] / divisor[divisor_start_idx];
			quotient[i] = factor;
			for (SizeT j = 0; j < m; ++j) {
				remainder[i+j] -= factor * divisor[divisor_start_idx+j];
			}
		}

		remainder.erase(remainder.begin(), remainder.end() - m + 1);
		return {quotient, remainder};
	}
}