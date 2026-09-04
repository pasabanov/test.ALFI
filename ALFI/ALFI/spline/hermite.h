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

		using Type = std::variant<typename Types::Classic,
								  typename Types::Cardinal,
								  typename Types::CatmullRom,
								  typename Types::KochanekBartels,
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

			/*
			* Identify the tangent method.
			*/
			enum class Method {
				Classic,
				Cardinal,
				CatmullRom,
				KochanekBartels,
				Akima,
				ModifiedAkima,
				Explicit
			};

			Method method = Method::Classic;

			Number cardinal_c;
			Number tension;
			Number bias;
			Number continuity;

			const Container<Number>* explicit_derivatives = nullptr;

			std::visit(util::misc::overload{
				[&](const typename Types::Classic&) { method = Method::Classic; },
				[&](const typename Types::Cardinal& c) { method = Method::Cardinal; cardinal_c = c.c; },
				[&](const typename Types::CatmullRom&) { method = Method::CatmullRom; },
				[&](const typename Types::KochanekBartels& kb) { method = Method::KochanekBartels; tension = kb.tension; bias = kb.bias; continuity = kb.continuity; },
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
							<< ": Explicit derivatives (of size "
							<< explicit_derivatives->size()
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
			Number clamped_left;
			Number clamped_right;

			std::visit(util::misc::overload{
				[&](const typename Boundaries::Linear&) { boundary_method = BoundaryMethod::Linear; },
				[&](const typename Boundaries::Quadratic&) { boundary_method = BoundaryMethod::Quadratic; polynomial_degree = 2; },
				[&](const typename Boundaries::Cubic&) { boundary_method = BoundaryMethod::Cubic; polynomial_degree = 3; },
				[&](const typename Boundaries::Polynomial& p) { boundary_method = BoundaryMethod::Polynomial; polynomial_degree = p.degree; },
				[&](const typename Boundaries::Clamped& c) { boundary_method = BoundaryMethod::Clamped; clamped_left = c.left; clamped_right = c.right; },
				[&](const typename Boundaries::Periodic&) { boundary_method = BoundaryMethod::Periodic; },
			}, boundaries_type);

			const bool periodic = boundary_method == BoundaryMethod::Periodic;

			/*
			* Periodic data is represented with the last point duplicated.
			* Therefore there are n - 1 unique points / intervals in one period.
			*/
			const SizeT period_points = periodic ? n - 1 : n;

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
			*
			* This is intentionally implemented directly instead of going through
			* simple_spline(), because here we only need P'(x_endpoint).
			*/
			const auto polynomial_endpoint_derivative =
				[&](bool right, SizeT degree) -> Number {
					if (degree == 0) {
						return 0;
					}

					const SizeT count = std::min<SizeT>(degree + 1, n);

					if (count < 2) {
						return 0;
					}

					const SizeT first = right ? n - count : 0;
					const SizeT last = first + count - 1;
					const SizeT r = right ? last : first;

					Number result = 0;

					/*
					* Lagrange basis derivative at x_r.
					*
					* For j != r:
					*
					* L'_j(x_r) =
					*   1 / (x_j - x_r)
					*   * product_{k != j,r}
					*     (x_r - x_k) / (x_j - x_k)
					*
					* Since sum_j L_j(x) = 1, the derivative of L_r
					* is minus the sum of the other basis derivatives.
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

			/*
			* Tangent arrays.
			*
			* m_in[i]  - tangent arriving at point i
			* m_out[i] - tangent leaving point i
			*
			* For every method except Kochanek-Bartels these are identical.
			*/
			Container<Number> m_in(n);
			Container<Number> m_out(n);

			/*
			* In periodic mode, coordinates of the wrapped neighbours have to be
			* shifted by +/- period. This matters when X is not uniformly spaced.
			*/
			const Number period = periodic ? X[n-1] - X[0] : 0;

			const auto wrapped_index = [&](std::ptrdiff_t i) -> SizeT {
				const std::ptrdiff_t m = static_cast<std::ptrdiff_t>(period_points);

				i %= m;
				if (i < 0) {
					i += m;
				}

				return static_cast<SizeT>(i);
			};

			const auto wrapped_x = [&](std::ptrdiff_t i) -> Number {
				if (!periodic) {
					return X[static_cast<SizeT>(i)];
				}

				const std::ptrdiff_t m = static_cast<std::ptrdiff_t>(period_points);

				const std::ptrdiff_t q =
					i >= 0
						? i / m
						: -((-i + m - 1) / m);

				const SizeT j = wrapped_index(i);

				return X[j] + q * period;
			};

			const auto wrapped_y = [&](std::ptrdiff_t i) -> Number {
				return Y[wrapped_index(i)];
			};

			/*
			* Periodic secant accessor.
			*/
			const auto periodic_delta = [&](std::ptrdiff_t i) -> Number {
				return delta[wrapped_index(i)];
			};

			/*
			* Akima's endpoint extension.
			*
			* Original Akima uses two extrapolated slopes on each side.
			* These correspond to quadratic extrapolation of the endpoint slopes.
			*
			* For very small data sets there are not enough distinct slopes,
			* so the available slope is simply continued.
			*/
			const auto akima_delta = [&](std::ptrdiff_t i) -> Number {
				if (periodic) {
					return periodic_delta(i);
				}
				if (i >= 0 && i < static_cast<std::ptrdiff_t>(interval_count)) {
					return delta[static_cast<SizeT>(i)];
				}
				if (interval_count == 1) {
					return delta[0];
				}
				if (i == -1) {
					return 2 * delta[0] - delta[1];
				}
				if (i == -2) {
					const Number d_minus_1 = 2 * delta[0] - delta[1];
					return 2 * d_minus_1 - delta[0];
				}
				if (i == static_cast<std::ptrdiff_t>(interval_count)) {
					return 2 * delta[interval_count-1] - delta[interval_count-2];
				}
				if (i == static_cast<std::ptrdiff_t>(interval_count) + 1) {
					const Number d_n = 2 * delta[interval_count-1] - delta[interval_count-2];
					return 2 * d_n - delta[interval_count-1];
				}
				/*
				* Should not be reached for the formulas below.
				*/
				return delta[std::clamp(static_cast<SizeT>(i), static_cast<SizeT>(0), interval_count - 1)];
			};

			const auto akima_tangent = [&](SizeT i, bool modified) -> Number {
				const Number d_im2 = akima_delta(static_cast<std::ptrdiff_t>(i) - 2);
				const Number d_im1 = akima_delta(static_cast<std::ptrdiff_t>(i) - 1);
				const Number d_i = akima_delta(static_cast<std::ptrdiff_t>(i));
				const Number d_ip1 = akima_delta(static_cast<std::ptrdiff_t>(i) + 1);

				Number w1 = std::abs(d_ip1 - d_i);
				Number w2 = std::abs(d_im1 - d_im2);

				if (modified) {
					w1 += std::abs(d_ip1 + d_i) / 2;
					w2 += std::abs(d_im1 + d_im2) / 2;
				}

				const Number w = w1 + w2;

				if (w == 0) {
					return (d_im1 + d_i) / 2;
				}

				return (w1 * d_im1 + w2 * d_i) / w;
			};

			/*
			* Compute tangents.
			*
			* In periodic mode we calculate the unique points 0 ... n-2 and
			* then copy the tangent(s) of point 0 to point n-1.
			*/
			if (method == Method::Explicit) {
				for (SizeT i = 0; i < n; ++i) {
					m_in[i] = (*explicit_derivatives)[i];
					m_out[i] = (*explicit_derivatives)[i];
				}
			} else {
				const SizeT tangent_count = periodic ? period_points : n;

				for (SizeT i = 0; i < tangent_count; ++i) {
					if (!periodic && (i == 0 || i == n - 1)) {
						continue;
					}

					switch (method) {
						case Method::Classic: {
							const Number left = periodic ? periodic_delta(static_cast<std::ptrdiff_t>(i) - 1) : delta[i-1];
							const Number right = periodic ? periodic_delta(static_cast<std::ptrdiff_t>(i)) : delta[i];
							const Number m = (left + right) / 2;
							m_in[i] = m;
							m_out[i] = m;
							break;
						}
						case Method::Cardinal: {
							const Number xp = wrapped_x(static_cast<std::ptrdiff_t>(i) - 1);
							const Number xn = wrapped_x(static_cast<std::ptrdiff_t>(i) + 1);
							const Number yp = wrapped_y(static_cast<std::ptrdiff_t>(i) - 1);
							const Number yn = wrapped_y(static_cast<std::ptrdiff_t>(i) + 1);
							const Number m = (1 - cardinal_c) * (yn - yp) / (xn - xp);
							m_in[i] = m;
							m_out[i] = m;
							break;
						}
						case Method::CatmullRom: {
							const Number xp = wrapped_x(static_cast<std::ptrdiff_t>(i) - 1);
							const Number xn = wrapped_x(static_cast<std::ptrdiff_t>(i) + 1);
							const Number yp = wrapped_y(static_cast<std::ptrdiff_t>(i) - 1);
							const Number yn = wrapped_y(static_cast<std::ptrdiff_t>(i) + 1);
							const Number m = (yn - yp) / (xn - xp);
							m_in[i] = m;
							m_out[i] = m;
							break;
						}
						case Method::KochanekBartels: {
							const Number d_left = periodic ? periodic_delta(static_cast<std::ptrdiff_t>(i) - 1) : delta[i-1];
							const Number d_right = periodic ? periodic_delta(static_cast<std::ptrdiff_t>(i)) : delta[i];

							const Number scale = (1 - tension) / 2;

							/*
							* Incoming tangent DS:
							*/
							m_in[i] =
								scale *
								(
									(1 + bias)
										* (1 - continuity)
										* d_left
									+
									(1 - bias)
										* (1 + continuity)
										* d_right
								);

							/*
							* Outgoing tangent DD:
							*/
							m_out[i] =
								scale *
								(
									(1 + bias)
										* (1 + continuity)
										* d_left
									+
									(1 - bias)
										* (1 - continuity)
										* d_right
								);

							break;
						}
						case Method::Akima:
							m_in[i] = akima_tangent(i, false);
							m_out[i] = m_in[i];
							break;
						case Method::ModifiedAkima:
							m_in[i] = akima_tangent(i, true);
							m_out[i] = m_in[i];
							break;
						case Method::Explicit:
							__builtin_unreachable();
					}
				}

				/*
				* Non-periodic endpoint conditions.
				*/
				if (!periodic) {
					switch (boundary_method) {
						case BoundaryMethod::Linear:
							m_out[0] = delta[0];
							m_in[n-1] = delta[n-2];
							break;
						case BoundaryMethod::Quadratic:
							m_out[0] = polynomial_endpoint_derivative(false, 2);
							m_in[n-1] = polynomial_endpoint_derivative(true, 2);
							break;
						case BoundaryMethod::Cubic:
							m_out[0] = polynomial_endpoint_derivative(false, 3);
							m_in[n-1] = polynomial_endpoint_derivative(true, 3);
							break;
						case BoundaryMethod::Polynomial:
							m_out[0] = polynomial_endpoint_derivative(false, polynomial_degree);
							m_in[n-1] = polynomial_endpoint_derivative(true, polynomial_degree);
							break;
						case BoundaryMethod::Clamped:
							m_out[0] = clamped_left;
							m_in[n-1] = clamped_right;
							break;
						case BoundaryMethod::Periodic:
							__builtin_unreachable();
					}

					/*
					* For ordinary Hermite variants the endpoint has only one
					* tangent. For KB, only the outgoing tangent at the left
					* endpoint and incoming tangent at the right endpoint are
					* relevant.
					*/
					if (method != Method::KochanekBartels) {
						m_in[0] = m_out[0];
						m_out[n-1] = m_in[n-1];
					}
				} else {
					/*
					* Duplicate the periodic endpoint.
					*/
					m_in[n-1] = m_in[0];
					m_out[n-1] = m_out[0];
				}
			}

			/*
			* Convert every Hermite segment to:
			*
			*   S(x) = a * z^3 + b * z^2 + c * z + d
			*
			* where z = x - X[i].
			*/
			Container<Number> coeffs;
			coeffs.resize(4 * interval_count);

			for (SizeT i = 0; i < interval_count; ++i) {
				const Number h = X[i+1] - X[i];

				const Number m0 = m_out[i];
				const Number m1 = m_in[i+1];

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