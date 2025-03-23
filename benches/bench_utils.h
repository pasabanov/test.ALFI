#pragma once

#include <random>

template <typename T>
T random_value(T min, T max) {
	static std::mt19937_64 generator(std::random_device{}());

	if constexpr (std::is_integral_v<T>) {
		static std::uniform_int_distribution<T> dist(min, max);
		return dist(generator);
	} else if constexpr (std::is_floating_point_v<T>) {
		static std::uniform_real_distribution<T> dist(min, max);
		return dist(generator);
	} else {
		static_assert(false, "Unsupported type: must be integral or floating-point");
		return {};
	}
}

template <typename F>
std::vector<double> apply_func(const std::vector<double>& input, F func) {
	std::vector<double> result(input.size());
	std::transform(input.begin(), input.end(), result.begin(), func);
	return result;
}