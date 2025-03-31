#include <iomanip>
#include <iostream>
#include <vector>

#include <ALFI/ratf.h>

constexpr double factorial(int n) {
	double result = 1;
	for (int i = 2; i <= n; ++i) {
		result *= i;
	}
	return result;
}

std::vector<double> maclaurin_exp(int order) {
	std::vector<double> coeffs(order + 1);
	for (int n = 0; n <= order; ++n) {
		coeffs[order-n] = 1.0 / factorial(n);
	}
	return coeffs;
}

std::vector<double> maclaurin_sin(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 1) {
			coeffs[order-n] = (n % 4 == 1 ? 1.0 : -1.0) / factorial(n);
		}
	}
	return coeffs;
}

std::vector<double> maclaurin_cos(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 0) {
			coeffs[order-n] = (n % 4 == 0 ? 1.0 : -1.0) / factorial(n);
		}
	}
	return coeffs;
}

std::vector<double> maclaurin_ln(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 1; n <= order; ++n) {
		coeffs[order-n] = (n % 2 == 1 ? 1.0 : -1.0) / n;
	}
	return coeffs;
}

std::vector<double> maclaurin_arctan(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 1) {
			coeffs[order-n] = (n % 4 == 1 ? 1.0 : -1.0) / n;
		}
	}
	return coeffs;
}

void print_vector(const std::vector<double>& poly, int precision = 6) {
	std::cout << std::fixed << std::setprecision(precision) << '[';
	for (size_t i = 0; i < poly.size(); ++i) {
		std::cout << poly[i] << (i + 1 < poly.size() ? ", " : "");
	}
	std::cout << ']';
}

void print_polynomial(const std::vector<double>& poly, int precision = 6) {
	std::cout << std::fixed << std::setprecision(precision);
	bool first = true;
	for (size_t i = 0; i < poly.size(); ++i) {
		const auto coeff = poly[i];
		if (coeff == 0) {
			continue;
		}

		if (!first) {
			std::cout << (coeff > 0 ? " + " : " - ");
		} else if (coeff < 0) {
			std::cout << '-';
		}

		const auto abs_coeff = std::abs(coeff);
		if (abs_coeff != 1 || i == poly.size() - 1) {
			std::cout << abs_coeff;
		}

		const auto power = poly.size() - 1 - i;
		if (power == 1) {
			std::cout << 'x';
		} else if (power > 1) {
			std::cout << "x^{" << power << '}';
		}

		first = false;
	}
}

void show_pade(const std::string& name, const std::vector<double>& series, size_t n, size_t m, int precision = 6) {
	const auto approx = alfi::ratf::pade(series, n, m);

	std::cout << "=== " << name << " (Padé [" << n << "/" << m << "]) ===\n";

	std::cout << "Maclaurin series:\n";
	std::cout << "  Coefficients: ";
	print_vector(series, precision);
	std::cout << "\n  Polynomial:   ";
	print_polynomial(series, precision);
	std::cout << "\n\n";

	std::cout << "Padé approximation:\n";
	std::cout << "  Numerator:    ";
	print_vector(approx.first, precision);
	std::cout << "\n  Denominator:  ";
	print_vector(approx.second, precision);
	std::cout << '\n';
	std::cout << "  Fraction:     \\frac{";
	print_polynomial(approx.first, precision);
	std::cout << "}{";
	print_polynomial(approx.second, precision);
	std::cout << "}\n\n";
}

int main() {
	const size_t order = 10;
	const size_t n = 4, m = 4;
	const int precision = 8;

	show_pade("e^x", maclaurin_exp(order), n, m, precision);
	show_pade("sin(x)", maclaurin_sin(order), n, m, precision);
	show_pade("cos(x)", maclaurin_cos(order), n, m, precision);
	show_pade("ln(1 + x)", maclaurin_ln(order), n, m, precision);
	show_pade("arctan(x)", maclaurin_arctan(order), n, m, precision);

	return 0;
}