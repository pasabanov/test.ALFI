#pragma once

#include <cmath>

#include "config.h"
#include "util/points.h"

namespace alfi::dist {
	enum class Type {
		GENERAL,
		UNIFORM,
		QUADRATIC,
		CUBIC,
		CHEBYSHEV,
		CHEBYSHEV_STRETCHED,
		CHEBYSHEV_AUGMENTED,
		CHEBYSHEV_2,
		CHEBYSHEV_3,
		CHEBYSHEV_3_STRETCHED,
		CHEBYSHEV_4,
		CHEBYSHEV_4_STRETCHED,
		CHEBYSHEV_ELLIPSE,
		CHEBYSHEV_ELLIPSE_STRETCHED,
		CHEBYSHEV_ELLIPSE_AUGMENTED,
		CHEBYSHEV_ELLIPSE_2,
		CHEBYSHEV_ELLIPSE_3,
		CHEBYSHEV_ELLIPSE_3_STRETCHED,
		CHEBYSHEV_ELLIPSE_4,
		CHEBYSHEV_ELLIPSE_4_STRETCHED,
		LOGISTIC,
		LOGISTIC_STRETCHED,
		ERF,
		ERF_STRETCHED,
	};

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> uniform(SizeT n, const Number& a, const Number& b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			points[i] = a + (b - a) * static_cast<Number>(i) / static_cast<Number>(n - 1);
		}
		return points;
	}

	/**
		@brief Generates a distribution of \p n points on the segment `[a, b]` using two quadratic polynomials.

		The following transform function \f(f\f):
		\f[
			f(x) = \begin{cases}
				(x+1)^2 - 1, & x \leq 0 \\
				-(x-1)^2 + 1, & x > 0
			\end{cases}
		\f]
		is applied to `n` points uniformly distributed on the segment `[-1, 1]`, mapping them onto the same segment.\n
		The resulting points are then linearly mapped to the target segment `[a, b]`.

		@param n number of points
		@param a left boundary of the segment
		@param b right boundary of the segment
		@return a container with \p n points distributed on the segment `[a, b]` according to the transform function
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> quadratic(SizeT n, const Number& a, const Number& b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 2 * static_cast<Number>(i) / (static_cast<Number>(n) - 1) - 1;
			const Number value = x <= 0 ? x * (x + 2) : -x * (x - 2);
			points[i] = a + (b - a) * (1 + value) / 2;
		}
		return points;
	}

	/**
		@brief Generates a distribution of \p n points on the segment `[a, b]` using cubic polynomial.

		The following transform function \f(f\f):
		\f[
			f(x) = -0.5x^3 + 1.5x
		\f]
		is applied to `n` points uniformly distributed on the segment `[-1, 1]`, mapping them onto the same segment.\n
		The resulting points are then linearly mapped to the target segment `[a, b]`.

		@param n number of points
		@param a left boundary of the segment
		@param b right boundary of the segment
		@return a container with \p n points distributed on the segment `[a, b]` according to the transform function
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> cubic(SizeT n, const Number& a, const Number& b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 2 * static_cast<Number>(i) / (static_cast<Number>(n) - 1) - 1;
			const Number value = (3 - x*x) * x / 2;
			points[i] = a + (b - a) * (1 + value) / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev(SizeT n, const Number& a, const Number& b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = std::cos(M_PI * (2 * (static_cast<Number>(n) - static_cast<Number>(i)) - 1) / (2 * static_cast<Number>(n)));
			points[i] = a + (x + 1) * (b - a) / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_stretched(SizeT n, const Number& a, const Number& b) {
		return points::stretched<Number,Container>(chebyshev(n, a, b), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_augmented(SizeT n, const Number& a, const Number& b) {
		if (n == 0) {
			return {};
		}
		if (n == 1) {
			return {(a+b)/2};
		}
		Container<Number> points(n);
		Container<Number> c = chebyshev(n - 2, a, b);
		std::move(c.begin(), c.end(), points.begin() + 1);
		points[0] = a;
		points[n-1] = b;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_2(SizeT n, const Number& a, const Number& b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 1 - std::cos(M_PI * static_cast<Number>(i) / (static_cast<Number>(n) - 1));
			points[i] = a + (b - a) * x / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_3(SizeT n, const Number& a, const Number& b) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 1 - std::cos(M_PI * static_cast<Number>(2*i) / static_cast<Number>(2*n - 1));
			points[i] = a + (b - a) * x / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_3_stretched(SizeT n, const Number& a, const Number& b) {
		return points::stretched<Number,Container>(chebyshev_3(n, a, b), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_4(SizeT n, const Number& a, const Number& b) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 1 - std::cos(M_PI * static_cast<Number>(2*i + 1) / static_cast<Number>(2*n - 1));
			points[i] = a + (b - a) * x / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_4_stretched(SizeT n, const Number& a, const Number& b) {
		return points::stretched<Number,Container>(chebyshev_4(n, a, b), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n / 2; ++i) {
			const Number theta = M_PI * (2 * static_cast<Number>(i) + 1) / (2 * static_cast<Number>(n));
			const Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = a + (1 - x) * (b - a) / 2;
			points[n-1-i] = a + (1 + x) * (b - a) / 2;
		}
		if (n % 2 == 1)
			points[n/2] = (a+b)/2;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_stretched(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		return points::stretched<Number,Container>(chebyshev_ellipse(n, a, b, ratio), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_augmented(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		if (n == 0) {
			return {};
		}
		if (n == 1) {
			return {(a+b)/2};
		}
		Container<Number> points(n);
		Container<Number> c = chebyshev_ellipse(n - 2, a, b, ratio);
		std::move(c.begin(), c.end(), points.begin() + 1);
		points[0] = a;
		points[n-1] = b;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_2(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n / 2; ++i) {
			const Number theta = M_PI * static_cast<Number>(i) / (static_cast<Number>(n) - 1);
			const Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = (1 - x) * (b - a) / 2 + a;
			points[n-1-i] = (1 + x) * (b - a) / 2 + a;
		}
		if (n % 2 == 1)
			points[n/2] = (a+b)/2;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_3(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number theta = M_PI * static_cast<Number>(2*i) / static_cast<Number>(2*n - 1);
			const Number x = (theta < M_PI/2 ? 1 : -1) / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = (1 - x) * (b - a) / 2 + a;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_3_stretched(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		return points::stretched<Number,Container>(chebyshev_ellipse_3(n, a, b, ratio), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_4(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number theta = M_PI * static_cast<Number>(2*i + 1) / static_cast<Number>(2*n - 1);
			const Number x = (theta < M_PI/2 ? 1 : -1) / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = (1 - x) * (b - a) / 2 + a;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_4_stretched(SizeT n, const Number& a, const Number& b, const Number& ratio) {
		return points::stretched<Number,Container>(chebyshev_ellipse_4(n, a, b, ratio), a, b);
	}

	/**
		@brief Generates a distribution of \p n points on the interval `(a, b)` using the logistic function.

		The following transform function \f(f\f):
		\f[
			f(x) = \frac{1}{1 + e^{-steepness \cdot x}}
		\f]
		is applied to `n` points uniformly distributed on the interval `(-1, 1)`, mapping them onto the `(0, 1)` interval.\n
		The resulting points are then linearly mapped to the target interval `(a, b)`.

		The slope of the transform function \f(f\f) at `x = 0` is determined by \p steepness and is given by:
		\f[
			\left. \dv{f}{x} \right\vert_{x=0} = \frac14 \cdot steepness
		\f]

		@note Extreme points don't lie on the interval `(a, b)` boundaries.

		@param n number of points
		@param a left boundary of the interval
		@param b right boundary of the interval
		@param steepness determines the slope of the transform function at `x = 0`
		@return a container with \p n points distributed on the interval `(a, b)` according to the transform function
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> logistic(SizeT n, const Number& a, const Number& b, const Number& steepness) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 2 * static_cast<Number>(i) / (static_cast<Number>(n) - 1) - 1;
			const Number logisticValue = 1 / (1 + std::exp(-steepness * x));
			points[i] = a + (b - a) * logisticValue;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> logistic_stretched(SizeT n, const Number& a, const Number& b, const Number& steepness) {
		if (n == 0)
			return {};
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		const Number stretch_factor = 1 - 2 / (1 + std::exp(steepness));
		for (SizeT i = 1; i < n - 1; ++i) {
			const Number x = 2 * static_cast<double>(i) / (static_cast<Number>(n) - 1) - 1;
			const Number logisticValue = 1 / (1 + std::exp(-steepness * x));
			const Number stretched = (logisticValue - 1 / (1 + std::exp(steepness))) / stretch_factor;
			points[i] = a + (b - a) * stretched;
		}
		points[0] = a;
		points[n-1] = b;
		return points;
	}

	/**
		@brief Generates a distribution of \p n points on the interval `(a, b)` using the error function.

		The following transform function \f(f\f):
		\f[
			f(x) = \operatorname{erf}(steepness \cdot x)
		\f]
		is applied to `n` points uniformly distributed on the interval `(-1, 1)`, mapping them onto the same interval.\n
		The resulting points are then linearly mapped to the target interval `(a, b)`.

		The slope of the transform function \f(f\f) at `x = 0` is determined by \p steepness and is given by:
		\f[
			\left. \dv{f}{x} \right\vert_{x=0} = \frac2\pi \cdot steepness
		\f]

		@note Extreme points don't lie on the interval `(a, b)` boundaries.

		@param n number of points
		@param a left boundary of the interval
		@param b right boundary of the interval
		@param steepness determines the slope of the transform function at `x = 0`
		@return a container with \p n points distributed on the interval `(a, b)` according to the transform function
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> erf(SizeT n, const Number& a, const Number& b, const Number& steepness) {
		if (n == 0)
			return {};
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 2 * static_cast<Number>(i) / (static_cast<Number>(n) - 1) - 1;
			const Number erf_value = std::erf(steepness * x);
			points[i] = a + (b - a) * (1 + erf_value) / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> erf_stretched(SizeT n, const Number& a, const Number& b, const Number& steepness) {
		return points::stretched<Number,Container>(erf(n, a, b, steepness), a, b);
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> of_type(Type type, SizeT n, const Number& a, const Number& b, const Number& parameter = NAN) {
		switch (type) {
		case Type::QUADRATIC:
			return quadratic(n, a, b);
		case Type::CUBIC:
			return cubic(n, a, b);
		case Type::CHEBYSHEV:
			return chebyshev(n, a, b);
		case Type::CHEBYSHEV_STRETCHED:
			return chebyshev_stretched(n, a, b);
		case Type::CHEBYSHEV_AUGMENTED:
			return chebyshev_augmented(n, a, b);
		case Type::CHEBYSHEV_2:
			return chebyshev_2(n, a, b);
		case Type::CHEBYSHEV_3:
			return chebyshev_3(n, a, b);
		case Type::CHEBYSHEV_3_STRETCHED:
			return chebyshev_3_stretched(n, a, b);
		case Type::CHEBYSHEV_4:
			return chebyshev_4(n, a, b);
		case Type::CHEBYSHEV_4_STRETCHED:
			return chebyshev_4_stretched(n, a, b);
		case Type::CHEBYSHEV_ELLIPSE:
			return chebyshev_ellipse(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_STRETCHED:
			return chebyshev_ellipse_stretched(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_AUGMENTED:
			return chebyshev_ellipse_augmented(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_2:
			return chebyshev_ellipse_2(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_3:
			return chebyshev_ellipse_3(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_3_STRETCHED:
			return chebyshev_ellipse_3_stretched(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_4:
			return chebyshev_ellipse_4(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_4_STRETCHED:
			return chebyshev_ellipse_4_stretched(n, a, b, parameter);
		case Type::LOGISTIC:
			return logistic(n, a, b, parameter);
		case Type::LOGISTIC_STRETCHED:
			return logistic_stretched(n, a, b, parameter);
		case Type::ERF:
			return erf(n, a, b, parameter);
		case Type::ERF_STRETCHED:
			return erf_stretched(n, a, b, parameter);
		case Type::UNIFORM:
		case Type::GENERAL:
		default:
			return uniform(n, a, b);
		}
	}
}