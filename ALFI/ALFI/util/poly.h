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
		auto p_start = std::find_if(p.begin(), p.end(), [&epsilon](Number v) { return std::abs(v) >= epsilon; });
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
			= std::find_if(divisor.begin(), divisor.end(), [&epsilon](Number v) { return std::abs(v) >= epsilon; });

		if (divisor_start == divisor.end() || dividend.size() < divisor.size()) {
			return {{}, dividend};
		}

		const auto divisor_start_idx = divisor_start - divisor.begin();

		Container<Number> quotient(dividend.size() - divisor.size() + 1, 0);
		Container<Number> remainder = dividend;

		for (SizeT i = 0; i <= dividend.size() - divisor.size(); ++i) {
			const Number factor = remainder[i] / divisor[divisor_start_idx];
			quotient[i] = factor;
			for (SizeT j = 0; j < divisor.size(); ++j) {
				remainder[i+j] -= factor * divisor[divisor_start_idx+j];
			}
		}

		remainder.erase(remainder.begin(), remainder.end() - divisor.size() + 1);
		return {quotient, remainder};
	}

	/**
		@brief Computes the extended Euclidean algorithm for polynomials.

		This function implements the extended Euclidean algorithm for polynomials, returning a tuple containing
		the greatest common divisor \f(r)\f (gcd) and the Bézout coefficients \f(s)\f and \f(t)\f such that:
		\f[
			a \cdot s + b \cdot t = r
		\f]
		where \f(a\f) and \f(b\f) are the input polynomials.

		The algorithm works with polynomials represented as containers of coefficients in descending order.
		The computation continues while the remainder \f(r)\f is non-empty and not equal (up to \p epsilon)
		to the zero polynomial (represented as a container with a single element close to zero).

		Note that the implementation avoids using separate normalization calls by checking that \f(r)\f is non-empty
		and not equivalent to zero in the loop condition.

		@param a the first polynomial (dividend)
		@param b the second polynomial (divisor)
		@param epsilon a threshold for treating coefficients as zero (default is machine epsilon)
		@return a tuple {r, s, t} where \f(r = \gcd(a, b)\f) and the Bézout identity \f(a \cdot s + b \cdot t = r\f) holds.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::tuple<Container<Number>,Container<Number>,Container<Number>>
	extended_euclid(const Container<Number>& a, const Container<Number>& b, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		Container<Number> old_r = a;
		Container<Number> r = b;

		Container<Number> old_s = {1};
		Container<Number> s = {0};
		Container<Number> old_t = {0};
		Container<Number> t = {1};

		while (!r.empty() && !(r.size() == 1 && std::abs(r[0]) < epsilon)) {
			const auto [q, new_r] = div(old_r, r, epsilon);

			old_r = std::move(r);
			r = std::move(new_r);

			const Container<Number> qs = mul(q, s);
			const auto new_s_size = std::max(old_s.size(), qs.size());
			Container<Number> new_s(new_s_size, 0);
			for (SizeT i = 0, offset = new_s_size - old_s.size(); i < old_s.size(); ++i) {
				new_s[offset+i] = old_s[i];
			}
			for (SizeT i = 0, offset = new_s_size - qs.size(); i < qs.size(); ++i) {
				new_s[offset+i] -= qs[i];
			}

			const Container<Number> qt = mul(q, t);
			const auto new_t_size = std::max(old_t.size(), qt.size());
			Container<Number> new_t(new_t_size, 0);
			for (SizeT i = 0, offset = new_t_size - old_t.size(); i < old_t.size(); ++i) {
				new_t[offset+i] = old_t[i];
			}
			for (SizeT i = 0, offset = new_t_size - qt.size(); i < qt.size(); ++i) {
				new_t[offset+i] -= qt[i];
			}

			old_s = std::move(s);
			s = std::move(new_s);
			old_t = std::move(t);
			t = std::move(new_t);
		}

		return std::make_tuple(old_r, old_s, old_t);
	}
}