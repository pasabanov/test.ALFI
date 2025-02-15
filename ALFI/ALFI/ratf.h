#pragma once

#include <iostream>
#include <cmath>

#include "config.h"

namespace alfi::ratf {
//	template <typename Number, template <typename, typename...> class Container = DefaultContainer>
//	using RationalFunction = std::pair<Container<Number>,Container<Number>>;

	template <typename Number, template <typename, typename...> class Container = DefaultContainer>
	class RationalFunction {
	public:
		Container<Number> num;
		Container<Number> den;

		RationalFunction() = default;

		RationalFunction(Container<Number> numerator, Container<Number> denominator) : num(std::move(numerator)), den(std::move(denominator)) {}

		RationalFunction(const RationalFunction& other) = default;
		RationalFunction(RationalFunction&& other) noexcept = default;

		// ReSharper disable once CppNonExplicitConvertingConstructor
		RationalFunction(std::pair<Container<Number>,Container<Number>> p) : num(std::move(p.first)), den(std::move(p.second)) {} // NOLINT(*-explicit-constructor)

		RationalFunction& operator=(const RationalFunction& other) = default;
		RationalFunction& operator=(RationalFunction&& other) noexcept = default;

		RationalFunction& operator=(std::pair<Container<Number>,Container<Number>> p) {
			num = std::move(p.first);
			den = std::move(p.second);
			return *this;
		}

		// ReSharper disable once CppNonExplicitConversionOperator
		operator std::pair<Container<Number>,Container<Number>>() const & { // NOLINT(*-explicit-constructor)
			return std::make_pair(num, den);
		}
		// ReSharper disable once CppNonExplicitConversionOperator
		operator std::pair<Container<Number>,Container<Number>>() && noexcept { // NOLINT(*-explicit-constructor)
			return std::make_pair(std::move(num), std::move(den));
		}
	};
}