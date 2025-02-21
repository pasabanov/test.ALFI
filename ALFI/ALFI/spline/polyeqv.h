#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>

#include "../config.h"
#include "../poly.h"
#include "../util/misc.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	class PolyEqvSpline {
	public:
		enum class PolynomialType {
			LAGRANGE,
			IMPROVED_LAGRANGE,
			NEWTON,

			Default = LAGRANGE,
		};

		enum class OptimizationType {
			ACCURACY,
			SPEED,

			Default = ACCURACY,
		};

		enum class EvaluationType {
			IGNORE_NANS_AND_PREVIOUS,
			IGNORE_NANS,
			NOT_IGNORE_NANS,

			Default = IGNORE_NANS_AND_PREVIOUS,
		};

		PolyEqvSpline() = default;
		PolyEqvSpline(const Container<Number>& X,
					  const Container<Number>& Y,
					  PolynomialType polynomial_type = PolynomialType::Default,
					  OptimizationType optimization_type = OptimizationType::Default) {
			construct(X, Y, polynomial_type, optimization_type);
		}

		PolyEqvSpline(const PolyEqvSpline& other) = default;
		PolyEqvSpline(PolyEqvSpline&& other) noexcept = default;

		PolyEqvSpline& operator=(const PolyEqvSpline& other) = default;
		PolyEqvSpline& operator=(PolyEqvSpline&& other) noexcept = default;

		void construct(const Container<Number>& X,
					   const Container<Number>& Y,
					   PolynomialType polynomial_type = PolynomialType::Default,
					   OptimizationType optimization_type = OptimizationType::Default,
					   EvaluationType evaluation_type = EvaluationType::Default) {
			_X = X;
			_polynomial_type = polynomial_type;
			_optimization_type = optimization_type;
			_evaluation_type = evaluation_type;

			const SizeT n = static_cast<SizeT>(X.size());

			if (n <= 3) {
				_coeffs = util::spline::simple_spline<Number,Container>(X, Y, n - 1);
				return;
			}

			_coeffs.resize(n * (n - 1)); // n coefficients for n-1 intervals

			switch (optimization_type) {
				case OptimizationType::ACCURACY: {
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
					for (SizeT i = 0; i < n - 1; ++i) {
						Container<Number> P;
						switch (polynomial_type) {
							case PolynomialType::LAGRANGE: P = poly::lagrange(X, Y); break;
							case PolynomialType::IMPROVED_LAGRANGE: P = poly::imp_lagrange(X, Y); break;
							case PolynomialType::NEWTON: P = poly::newton(X, Y); break;
						}
						std::move(P.begin(), P.end(), _coeffs.begin() + (i * n));
						_coeffs[i*n+n-1] = Y[i];
					}
					break;
				}
				case OptimizationType::SPEED: {
					Container<Number> P;
					switch (polynomial_type) {
						case PolynomialType::LAGRANGE: P = poly::lagrange(X, Y); break;
						case PolynomialType::IMPROVED_LAGRANGE: P = poly::imp_lagrange(X, Y); break;
						case PolynomialType::NEWTON: P = poly::newton(X, Y); break;
					}

					const static auto binomials = [](SizeT m) {
						Container<Container<Number>> C(m + 1);
						for (SizeT i = 0; i <= m; ++i) {
							C[i].resize(i + 1);
							C[i][0] = C[i][i] = 1;
							for (SizeT j = 1; j <= i / 2; ++j) {
								C[i][j] = C[i][i - j] = C[i - 1][j - 1] + C[i - 1][j];
							}
						}
						return C;
					};

					const auto C = binomials(n - 1);

#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
					for (SizeT i = 0; i < n - 1; ++i) {
						_coeffs[i*n] = P[0];
						for (SizeT k = 1; k < n - 1; ++k) {
							Number r = 0;
							for (SizeT j = 0; j < k; ++j) {
								r += (((k - j) & 1 == 0) ? 1 : -1) * _coeffs[i*n+j] * C[n-1-j][k-j];
								r *= X[i];
							}
							_coeffs[i*n+k] = P[k] - r;
						}
						_coeffs[i*n+n-1] = Y[i];
					}
					break;
				}
			}
		}

		Number eval(Number const x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, SizeT segment) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}

			segment = std::max(std::min(segment, _X.size() - 2), 0);

			x = x - _X[segment];

			const SizeT n = _X.size();

			Number result = _coeffs[segment*n];

			if (std::isnan(result)) {
				switch (_evaluation_type) {
					case EvaluationType::IGNORE_NANS_AND_PREVIOUS:
					case EvaluationType::IGNORE_NANS:
						result = 0;
						break;
					case EvaluationType::NOT_IGNORE_NANS:
						return result;
				}
			}

			for (SizeT i = 1; i < n; ++i) {
				result *= x;
				Number current = _coeffs[segment*n+i];
				if (std::isnan(current)) {
					switch (_evaluation_type) {
						case EvaluationType::IGNORE_NANS_AND_PREVIOUS:
							result = 0;
							current = 0;
							break;
						case EvaluationType::IGNORE_NANS:
							current = 0;
							break;
						case EvaluationType::NOT_IGNORE_NANS:
							return current;
					}
				}
				result += current;
			}

			return result;
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container<Number> result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number evalx = xx[i];
					while (i_x + 1 < _X.size() && evalx >= _X[i_x+1])
						++i_x;
					result[i] = eval(evalx, i_x);
				}
			} else {
				for (SizeT i = 0; i < xx.size(); ++i) {
					result[i] = eval(xx[i]);
				}
			}
			return result;
		}

		Number operator()(Number x) const {
			return eval(x);
		}
		Container<Number> operator()(const Container<Number>& xx) const {
			return eval(xx);
		}

		PolynomialType polynomial_type() const {
			return _polynomial_type;
		}

		OptimizationType optimization_type() const {
			return _optimization_type;
		}

		const Container<Number>& X() const {
			return _X;
		}
		Container<Number>&& X() && {
			return std::move(_X);
		}

		const Container<Number>& coeffs() const {
			return _coeffs;
		}
		Container<Number>&& coeffs() && {
			return std::move(_coeffs);
		}

		static std::pair<SizeT, SizeT> segment(SizeT index) {
			return {_X.size()*index, _X.size()*(index+1)};
		}

	private:
		PolynomialType _polynomial_type = PolynomialType::Default;
		OptimizationType _optimization_type = OptimizationType::Default;
		EvaluationType _evaluation_type = EvaluationType::Default;
		Container<Number> _X {};
		Container<Number> _coeffs {};
	};
}