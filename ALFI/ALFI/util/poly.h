#pragma once

#include <algorithm>
#include <limits>

#include "../config.h"

namespace alfi::util::poly {
	/**
		@brief Normalizes a polynomial by removing leading coefficients close to zero.

		This function trims coefficients from the beginning of the polynomial that are smaller than a given threshold \p epsilon.
		If the polynomial becomes empty after trimming, it adds to the container a single zero.

		@param p the polynomial to normalize (a container of coefficients)
		@param epsilon the threshold to treat coefficients as zero (default is machine epsilon)
		@return the normalized polynomial with leading zero coefficients removed
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

		Given two polynomials \f(p_1\f) and \f(p_2\f), represented as containers of coefficients in descending order,
		this function computes their product using the standard convolution formula:
		\f[
			(p_1 \cdot p_2)[k] = \sum_{i+j=k} p_1[i] p_2[j]
		\f]
		where indices correspond to powers of the variable.

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

		This function divides the polynomial \f(D\f) (dividend) by the polynomial \f(d\f) (divisor),
		returning the quotient \f(Q\f) and remainder \f(R\f) such that:
		\f[
			D(x) = Q(x) \cdot d(x) + R(x)
		\f]
		with \f(deg(R) < deg(d)\f) or \f(R = 0\f).

		Polynomials are represented as containers of coefficients in descending order.
		If the divisor is the zero polynomial, the function returns an empty quotient
		and the original dividend as the remainder.

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

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::tuple<Container<Number>,Container<Number>,Container<Number>> extended_euclid(const Container<Number>& a, const Container<Number>& b, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		// Копии входных многочленов
		Container<Number> old_r = a;
		Container<Number> r = b;

		// Начальные коэффициенты Bézout'а
		Container<Number> old_s = {1};
		Container<Number> s = {0};
		Container<Number> old_t = {0};
		Container<Number> t = {1};

		// Нормализуем входные многочлены
		// normalize(old_r, epsilon);
		// normalize(r, epsilon);

		// Основной цикл алгоритма
		// Если r равен нулевому многочлену (представленному как {0}), то завершаем цикл
		while (!r.empty() && !(r.size() == 1 && std::abs(r[0]) < epsilon)) {
			// Делим old_r на r: old_r = q * r + remainder
			auto [q, rem] = div(old_r, r, epsilon);

			// Обновляем (old_r, r): (old_r, r) := (r, rem)
			old_r = std::move(r);
			r = std::move(rem);
			// Если r оказался пустым, нормализуем его до {0}
			// normalize(r, epsilon);

			// Вычисляем новое s: new_s = old_s - q * s
			Container<Number> qs = mul(q, s);
			auto new_size = std::max(old_s.size(), qs.size());
			Container<Number> new_s(new_size, 0);
			auto offset = new_size - old_s.size();
			for (SizeT i = 0; i < old_s.size(); ++i) {
				new_s[offset + i] = old_s[i];
			}
			offset = new_size - qs.size();
			for (SizeT i = 0; i < qs.size(); ++i) {
				new_s[offset + i] -= qs[i];
			}

			// Аналогично для t: new_t = old_t - q * t
			Container<Number> qt = mul(q, t);
			new_size = std::max(old_t.size(), qt.size());
			Container<Number> new_t(new_size, 0);
			offset = new_size - old_t.size();
			for (SizeT i = 0; i < old_t.size(); ++i) {
				new_t[offset + i] = old_t[i];
			}
			offset = new_size - qt.size();
			for (SizeT i = 0; i < qt.size(); ++i) {
				new_t[offset + i] -= qt[i];
			}

			// Обновляем s и t: (old_s, s) := (s, new_s) и (old_t, t) := (t, new_t)
			old_s = std::move(s);
			s = std::move(new_s);
			old_t = std::move(t);
			t = std::move(new_t);

			// Нормализуем s и t
			// normalize(s, epsilon);
			// normalize(t, epsilon);
		}

		// По завершении алгоритма old_r содержит gcd, а old_s и old_t – коэффициенты Bézout'а.
		// normalize(old_r, epsilon);
		// normalize(old_s, epsilon);
		// normalize(old_t, epsilon);
		return std::make_tuple(old_r, old_s, old_t);
	}
}