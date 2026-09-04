#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>
#include <variant>

#include "../config.h"
#include "../util/misc.h"
#include "../util/spline.h"

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

		using Type = std::variant<typename Types::Classic,
								  typename Types::Cardinal,
								  typename Types::CatmullRom,
								  // typename Types::Pchip,
								  // typename Types::Hyman,
								  // typename Types::Steffen,
								  typename Types::Akima,
								  typename Types::ModifiedAkima,
								  typename Types::Explicit>;

		struct Boundaries final {
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

		using BoundariesType = std::variant<typename Boundaries::Linear,
											typename Boundaries::Quadratic,
											typename Boundaries::Cubic,
											typename Boundaries::Polynomial,
											typename Boundaries::Clamped,
											typename Boundaries::Periodic>;

		static Container<Number> compute_coeffs(
				const Container<Number>& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Returning an empty array..." << std::endl;
				return {};
			}

			const auto n = X.size();

			if (n <= 1) {
				return util::spline::simple_spline<Number,Container>(X, Y, 3);
			}

			enum class Method {
				Classic,
				Cardinal,
				CatmullRom,
				Akima,
				ModifiedAkima,
				Explicit
			};

			Method method = Method::Classic;

			Number cardinal_c {};

			const Container<Number>* explicit_derivatives = nullptr;

			std::visit(util::misc::overload{
				[&](const typename Types::Classic&) { method = Method::Classic; },
				[&](const typename Types::Cardinal& c) { method = Method::Cardinal; cardinal_c = c.c; },
				[&](const typename Types::CatmullRom&) { method = Method::CatmullRom; },
				[&](const typename Types::Akima&) { method = Method::Akima; },
				[&](const typename Types::ModifiedAkima&) { method = Method::ModifiedAkima; },
				[&](const typename Types::Explicit& e) { method = Method::Explicit; explicit_derivatives = &e.derivatives; },
			}, type);

			/*
			* Explicit derivatives completely override boundary handling.
			*/
			if (method == Method::Explicit) {
				if (explicit_derivatives->size() != n) {
					std::cerr << "Error in function " << __FUNCTION__
							<< ": Explicit derivatives (of size " << explicit_derivatives->size()
							<< ") and points (of size " << n
							<< ") are not the same size. Returning an empty array..."
							<< std::endl;
					return {};
				}
			}

			enum class BoundaryMethod {
				Linear,
				Quadratic,
				Cubic,
				Polynomial,
				Clamped,
				Periodic
			};

			BoundaryMethod boundary_method = BoundaryMethod::Linear;
			SizeT polynomial_degree = 1;
			Number clamped_left {};
			Number clamped_right {};

			std::visit(util::misc::overload{
				[&](const typename Boundaries::Linear&) { boundary_method = BoundaryMethod::Linear; },
				[&](const typename Boundaries::Quadratic&) { boundary_method = BoundaryMethod::Quadratic; polynomial_degree = 2; },
				[&](const typename Boundaries::Cubic&) { boundary_method = BoundaryMethod::Cubic; polynomial_degree = 3; },
				[&](const typename Boundaries::Polynomial& p) { boundary_method = BoundaryMethod::Polynomial; polynomial_degree = p.degree; },
				[&](const typename Boundaries::Clamped& c) { boundary_method = BoundaryMethod::Clamped; clamped_left = c.left; clamped_right = c.right; },
				[&](const typename Boundaries::Periodic&) { boundary_method = BoundaryMethod::Periodic; },
			}, boundaries_type);

			const SizeT interval_count = n - 1;

			/*
			* Secant slopes.
			*
			* In periodic mode d[n-2] is the slope of the closing interval,
			* which is already represented by X[n-1] and Y[n-1].
			*/
			Container<Number> delta(interval_count);

			for (SizeT i = 0; i < interval_count; ++i) {
				delta[i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
			}

			/*
			* Endpoint derivative obtained by differentiating an interpolating
			* polynomial through the first/last degree+1 points.
			*/
			const auto polynomial_endpoint_derivative =
				[&](bool left, SizeT degree) -> Number {
					if (degree == 0) {
						return 0;
					}

					const SizeT count = std::min<SizeT>(degree + 1, n);

					if (count < 2) {
						return 0;
					}

					const SizeT first = left ? 0 : n - count;
					const SizeT last = first + count - 1;
					const SizeT r = left ? first : last;

					Number result = 0;

					/*
					* Lagrange basis derivative at x_r.
					*
					* L'_j(x_r) = 1 / (x_j - x_r) * product_{k != j,r}{(x_r - x_k) / (x_j - x_k)}
					*/
					for (SizeT j = first; j <= last; ++j) {
						if (j == r) {
							continue;
						}
						Number basis_derivative = 1 / (X[j] - X[r]);
						for (SizeT k = first; k <= last; ++k) {
							if (k == j || k == r) {
								continue;
							}
							basis_derivative *= (X[r] - X[k]) / (X[j] - X[k]);
						}
						result += (Y[j] - Y[r]) * basis_derivative;
					}
					return result;
				};

			Container<Number> derivatives(n);

			/*
			* Compute tangents.
			*
			* In periodic mode we calculate the unique points 0 ... n-2 and
			* then copy the tangent(s) of point 0 to point n-1.
			*/
			if (method == Method::Explicit) {
				derivatives = *explicit_derivatives;
			} else {
				const SizeT end = boundary_method == BoundaryMethod::Periodic ? n : n - 1;
				for (SizeT i = boundary_method == BoundaryMethod::Periodic ? 0 : 1; i < end; ++i) {
					const SizeT prev_i = i == 0 ? n - 1 : i - 1;
					const SizeT next_i = i + 1 == n ? 0 : i + 1;

					switch (method) {
						case Method::Classic: {
							derivatives[i] = (delta[prev_i] + delta[i]) / 2;
							break;
						}
						case Method::Cardinal: {
							derivatives[i] = (1 - cardinal_c) * (Y[next_i] - Y[prev_i]) / (X[next_i] - X[prev_i]);
							break;
						}
						case Method::CatmullRom: {
							derivatives[i] = (Y[next_i] - Y[prev_i]) / (X[next_i] - X[prev_i]);
							break;
						}
						case Method::Akima:
						case Method::ModifiedAkima: {
							/*
							* Akima's endpoint extension.
							*
							* Original Akima uses two extrapolated slopes on each side.
							* These correspond to quadratic extrapolation of the endpoint slopes.
							*
							* Periodic Akima uses nodes from the opposite end of the array
							* as neighboring nodes, so no extrapolation is required.
							*/
							Number d_im2;
							Number d_im1;
							Number d_i;
							Number d_ip1;

							if (i == 0) {
								d_im2 = 3 * delta[0] - 2 * delta[1];
								d_im1 = 2 * delta[0] - delta[1];
								d_i = delta[0];
								d_ip1 = delta[1];
							} else if (i == 1) {
								d_im2 = delta[0];
								d_im1 = delta[1];
								d_i = delta[1];
								d_ip1 = delta[2];
							} else if (i == n - 2) {
								d_im2 = delta[i-2];
								d_im1 = delta[i-1];
								d_i = delta[i];
								d_ip1 = 2 * delta[i] - delta[i-1];
							} else if (i == n - 1) {
								d_im2 = delta[i-2];
								d_im1 = delta[i-1];
								d_i = 2 * delta[i-1] - delta[i-2];
								d_ip1 = 3 * delta[i-1] - 2 * delta[i-2];
							} else {
								d_im2 = delta[i-2];
								d_im1 = delta[i-1];
								d_i = delta[i];
								d_ip1 = delta[i+1];
							}

							Number w1 = std::abs(d_ip1 - d_i);
							Number w2 = std::abs(d_im1 - d_im2);

							if (method == Method::ModifiedAkima) {
								w1 += std::abs(d_ip1 + d_i) / 2;
								w2 += std::abs(d_im1 + d_im2) / 2;
							}

							const Number w = w1 + w2;

							if (w == 0) {
								derivatives[i] = (d_im1 + d_i) / 2;
							} else {
								derivatives[i] = (w1 * d_im1 + w2 * d_i) / w;
							}
							break;
						}
						case Method::Explicit:
							__builtin_unreachable();
					}
				}

				switch (boundary_method) {
					case BoundaryMethod::Linear:
						derivatives[0] = delta[0];
						derivatives[n-1] = delta[n-2];
						break;
					case BoundaryMethod::Quadratic:
						derivatives[0] = polynomial_endpoint_derivative(true, 2);
						derivatives[n-1] = polynomial_endpoint_derivative(false, 2);
						break;
					case BoundaryMethod::Cubic:
						derivatives[0] = polynomial_endpoint_derivative(true, 3);
						derivatives[n-1] = polynomial_endpoint_derivative(false, 3);
						break;
					case BoundaryMethod::Polynomial:
						derivatives[0] = polynomial_endpoint_derivative(true, polynomial_degree);
						derivatives[n-1] = polynomial_endpoint_derivative(false, polynomial_degree);
						break;
					case BoundaryMethod::Clamped:
						derivatives[0] = clamped_left;
						derivatives[n-1] = clamped_right;
						break;
					case BoundaryMethod::Periodic:
						// do nothing
						break;
				}
			}

			/*
			* Convert every Hermite segment to:
			*
			* S(x) = a * z^3 + b * z^2 + c * z + d
			*
			* where z = x - X[i].
			*/
			Container<Number> coeffs;
			coeffs.resize(4 * interval_count);

			for (SizeT i = 0; i < interval_count; ++i) {
				const Number h = X[i+1] - X[i];

				const Number m0 = derivatives[i];
				const Number m1 = derivatives[i+1];

				const Number y0 = Y[i];
				const Number y1 = Y[i+1];

				const Number h2 = h * h;
				const Number h3 = h2 * h;

				const Number a = (2 * y0 - 2 * y1 + h * (m0 + m1)) / h3;

				const Number b = (-3 * y0 + 3 * y1 - h * (2 * m0 + m1)) / h2;

				const Number c = m0;
				const Number d = y0;

				coeffs[4*i+0] = a;
				coeffs[4*i+1] = b;
				coeffs[4*i+2] = c;
				coeffs[4*i+3] = d;
			}

			return coeffs;
		}

		HermiteSpline() = default;

		template <typename ContainerXType>
		HermiteSpline(
				ContainerXType&& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			construct(std::forward<ContainerXType>(X), Y, type, boundaries_type);
		}

		HermiteSpline(const HermiteSpline& other) = default;
		HermiteSpline(HermiteSpline&& other) noexcept = default;

		HermiteSpline& operator=(const HermiteSpline& other) = default;
		HermiteSpline& operator=(HermiteSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(
				ContainerXType&& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}
			auto coeffs = compute_coeffs(X, Y, type, boundaries_type);
			if (coeffs.empty() && !X.empty()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": failed to construct coefficients. Not changing the object..." << std::endl;
				return;
			}
			_X = std::forward<ContainerXType>(X);
			_coeffs = std::move(coeffs);
		}

		Number eval(const Number& x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(const Number& x, SizeT segment) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}
			segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			const Number x_seg = x - _X[segment];
			return ((_coeffs[4*segment] * x_seg + _coeffs[4*segment+1]) * x_seg + _coeffs[4*segment+2]) * x_seg + _coeffs[4*segment+3];
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