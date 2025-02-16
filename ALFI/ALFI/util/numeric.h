#pragma once

#include <cmath>

#include "../config.h"

namespace alfi::util::numeric {
	template <typename Number = DefaultNumber>
	bool are_equal(Number a, Number b, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		return std::abs(a - b) <= epsilon || std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * epsilon;
	}

	template <typename Number, typename ExponentType>
	Number binpow(Number x, ExponentType n) {
		if constexpr (std::is_signed_v<ExponentType>) {
			if (n < 0) {
				return 1 / binpow(x, -n);
			}
		}
		Number result = 1;
		while (n > 0) {
			if ((n & 1) == 1) {
				result *= x;
			}
			x *= x;
			n >>= 1;
		}
		return result;
	}
}