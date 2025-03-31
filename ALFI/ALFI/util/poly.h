#pragma once

#include <algorithm>
#include <limits>

#include "../config.h"

namespace alfi::util::poly {
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

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::pair<Container<Number>,Container<Number>> div(const Container<Number>& dividend, const Container<Number>& divisor, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		const auto divisor_start = std::find_if(divisor.begin(), divisor.end(), [&epsilon](Number v) { return std::abs(v) > epsilon; });

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

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::tuple<Container<Number>,Container<Number>,Container<Number>> extended_euclid(
			Container<Number> a,
			Container<Number> b,
			Number epsilon = std::numeric_limits<Number>::epsilon(),
			SizeT min_r_degree = 0
	) {
		Container<Number> old_r = std::move(a);
		Container<Number> r = std::move(b);

		Container<Number> old_s = {1};
		Container<Number> s = {0};
		Container<Number> old_t = {0};
		Container<Number> t = {1};

		normalize(old_r);
		normalize(r);

		if (old_r.size() < r.size()) {
			std::swap(old_r, r);
			std::swap(old_s, s);
			std::swap(old_t, t);
		}

		while (!(r.size() == 1 && std::abs(r[0]) <= epsilon) && old_r.size() > min_r_degree + 1) {
			auto [q, new_r] = div(old_r, r, epsilon);

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

			old_r = std::move(r);
			r = std::move(new_r);
			old_s = std::move(s);
			s = std::move(new_s);
			old_t = std::move(t);
			t = std::move(new_t);

			normalize(old_r);
			normalize(r);
			normalize(old_s);
			normalize(s);
			normalize(old_t);
			normalize(t);
		}

		return std::make_tuple(old_r, old_s, old_t);
	}
}