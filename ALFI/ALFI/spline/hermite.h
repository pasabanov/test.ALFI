#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util/misc.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	class HermiteSpline {
	public:
		struct Types final {
			struct Classic final {};
			struct Cardinal final {
				explicit Cardinal(Number c) : c(std::move(c)) {}
				Number c;
			};
			struct CatmullRom final {};
			struct KochanekBartels final {
				explicit KochanekBartels(Number tension, Number bias, Number continuity) : tension(std::move(tension)), bias(std::move(bias)), continuity(std::move(continuity)) {}
				Number tension, bias, continuity;
			};
			// struct Pchip final {}; // TODO ?
			// struct Hyman final {}; // TODO ?
			// struct Steffen final {}; // TODO ?
			struct Akima final {};
			struct ModifiedAkima final {};
			struct Explicit final {
				explicit Explicit(Container<Number> derivatives) : derivatives(std::move(derivatives)) {}
				Container<Number> derivatives;
			};
			using Default = Classic;
		};

		struct Boundary final {
			struct Linear final {};
			struct Quadratic final {};
			struct Cubic final {};
			struct Polynomial final {
				explicit Polynomial(SizeT degree) : degree(std::move(degree)) {}
				SizeT degree;
			};
			struct Clamped final {
				explicit Clamped(Number left, Number right) : left(std::move(left)), right(std::move(right)) {}
				Number left, right;
			};
			struct Periodic final {};
			using Default = Linear;
		};

		static Container<Number> compute_coeffs(const Container<Number>& X, const Container<Number>& Y) {
			// if (X.size() != Y.size()) {
			// 	std::cerr << "Error in function " << __FUNCTION__
			// 			  << ": Vectors X (of size " << X.size()
			// 			  << ") and Y (of size " << Y.size()
			// 			  << ") are not the same size. Returning an empty array..." << std::endl;
			// 	return {};
			// }
			//
			// const auto n = X.size();
			//
			// if (n <= 2) {
			// 	if (n == 0) {
			// 		return {};
			// 	} else if (n == 1) {
			// 		return {Y[0]};
			// 	} else {
			// 		return {(Y[1] - Y[0]) / (X[1] - X[0]), Y[0]};
			// 	}
			// }
			//
			// Container<Number> coeffs;
			//
			// coeffs.resize(2 * (n - 1));
			//
			// for (SizeT i = 0; i < n - 1; ++i) {
			// 	coeffs[2*i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
			// 	coeffs[2*i+1] = Y[i];
			// }
			//
			// return coeffs;
		}

		HermiteSpline() = default;

		template <typename ContainerXType>
		HermiteSpline(ContainerXType&& X, const Container<Number>& Y) {
			construct(std::forward<ContainerXType>(X), Y);
		}

		HermiteSpline(const HermiteSpline& other) = default;
		HermiteSpline(HermiteSpline&& other) noexcept = default;

		HermiteSpline& operator=(const HermiteSpline& other) = default;
		HermiteSpline& operator=(HermiteSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const Container<Number>& Y) {
			// if (X.size() != Y.size()) {
			// 	std::cerr << "Error in function " << __FUNCTION__
			// 			  << ": Vectors X (of size " << X.size()
			// 			  << ") and Y (of size " << Y.size()
			// 			  << ") are not the same size. Doing nothing..." << std::endl;
			// 	return;
			// }
			// auto coeffs = compute_coeffs(X, Y);
			// if (coeffs.empty() && !X.empty()) {
			// 	std::cerr << "Error in function " << __FUNCTION__
			// 			  << ": failed to construct coefficients. Not changing the object..." << std::endl;
			// 	return;
			// }
			// _X = std::forward<ContainerXType>(X);
			// _coeffs = std::move(coeffs);
		}

		Number eval(const Number& x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(const Number& x, SizeT segment) const {
			// if (_coeffs.empty()) {
			// 	return NAN;
			// } else if (_coeffs.size() == 1) {
			// 	return _coeffs[0];
			// }
			// segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			// const Number x_seg = x - _X[segment];
			// return _coeffs[2*segment] * x_seg + _coeffs[2*segment+1];
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container<Number> result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number& x = xx[i];
					while (i_x + 1 < _X.size() && x >= _X[i_x+1])
						++i_x;
					result[i] = eval(x, i_x);
				}
			} else {
				for (SizeT i = 0; i < xx.size(); ++i) {
					result[i] = eval(xx[i]);
				}
			}
			return result;
		}

		Number operator()(const Number& x) const {
			return eval(x);
		}
		Container<Number> operator()(const Container<Number>& xx) const {
			return eval(xx);
		}

		const Container<Number>& X() const & {
			return _X;
		}
		Container<Number>&& X() && {
			return std::move(_X);
		}

		const Container<Number>& coeffs() const & {
			return _coeffs;
		}
		Container<Number>&& coeffs() && {
			return std::move(_coeffs);
		}

		static std::pair<SizeT, SizeT> segment(SizeT index) {
			return {2*index, 2*(index+1)};
		}

	private:
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}