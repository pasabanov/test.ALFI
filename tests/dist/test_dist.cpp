#include <ALFI/dist.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/dist/dist.toml";

const auto test_data = toml::parse_file(test_data_path);

const auto& mapping_intervals = test_data["mapping_intervals"].ref<toml::array>();

void test_distribution(const char*const name, const alfi::dist::Type type, double epsilon) {
	const auto& test_cases = test_data[name]["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto n = test_case["n"].ref<int64_t>();
		const auto a = test_case["a"].value<double>().value();
		const auto b = test_case["b"].value<double>().value();

		const auto parameter = test_case["ratio"].value_or(test_case["steepness"].value_or(static_cast<double>(NAN)));

		const auto expected = to_vector<double>(test_case["expected"].ref<toml::array>());
		const auto generated = of_type(type, n, a, b, parameter);
		expect_eq(generated, expected, epsilon);

		mapping_intervals.for_each([&](const toml::array& interval) {
			const auto c = interval[0].value<double>().value();
			const auto d = interval[1].value<double>().value();
			const auto mapped_expected = alfi::points::lin_mapped(expected, a, b, c, d);
			const auto mapped_generated = of_type(type, n, c, d, parameter);
			expect_eq(mapped_generated, mapped_expected, epsilon);
		});
	});
}

TEST(DistributionsTest, Uniform) {
	test_distribution("uniform", alfi::dist::Type::UNIFORM, 1e-15);
}

TEST(DistributionsTest, Quadratic) {
	test_distribution("quadratic", alfi::dist::Type::QUADRATIC, 1e-15);
}

TEST(DistributionsTest, Cubic) {
	test_distribution("cubic", alfi::dist::Type::CUBIC, 1e-15);
}

TEST(DistributionsTest, Chebyshev) {
	test_distribution("chebyshev", alfi::dist::Type::CHEBYSHEV, 1e-15);
}

TEST(DistributionsTest, ChebyshevStretched) {
	test_distribution("chebyshev_stretched", alfi::dist::Type::CHEBYSHEV_STRETCHED, 1e-11);
}

TEST(DistributionsTest, ChebyshevAugmented) {
	test_distribution("chebyshev_augmented", alfi::dist::Type::CHEBYSHEV_AUGMENTED, 1e-15);
}

TEST(DistributionsTest, Chebyshev2) {
	test_distribution("chebyshev_2", alfi::dist::Type::CHEBYSHEV_2, 1e-11);
}

TEST(DistributionsTest, Chebyshev3) {
	test_distribution("chebyshev_3", alfi::dist::Type::CHEBYSHEV_3, 1e-15);
}

TEST(DistributionsTest, Chebyshev3Stretched) {
	test_distribution("chebyshev_3_stretched", alfi::dist::Type::CHEBYSHEV_3_STRETCHED, 1e-15);
}

TEST(DistributionsTest, Chebyshev4) {
	test_distribution("chebyshev_4", alfi::dist::Type::CHEBYSHEV_4, 1e-15);
}

TEST(DistributionsTest, Chebyshev4Stretched) {
	test_distribution("chebyshev_4_stretched", alfi::dist::Type::CHEBYSHEV_4_STRETCHED, 1e-15);
}

TEST(DistributionsTest, ChebyshevEllipse) {
	test_distribution("chebyshev_ellipse", alfi::dist::Type::CHEBYSHEV_ELLIPSE, 1e-14);
}

TEST(DistributionsTest, ChebyshevEllipseStretched) {
	test_distribution("chebyshev_ellipse_stretched", alfi::dist::Type::CHEBYSHEV_ELLIPSE_STRETCHED, 1e-11);
}

TEST(DistributionsTest, ChebyshevEllipseAugmented) {
	test_distribution("chebyshev_ellipse_augmented", alfi::dist::Type::CHEBYSHEV_ELLIPSE_AUGMENTED, 1e-14);
}

TEST(DistributionsTest, ChebyshevEllipse2) {
	test_distribution("chebyshev_ellipse_2", alfi::dist::Type::CHEBYSHEV_ELLIPSE_2, 1e-14);
}

TEST(DistributionsTest, ChebyshevEllipse3) {
	test_distribution("chebyshev_ellipse_3", alfi::dist::Type::CHEBYSHEV_ELLIPSE_3, 1e-13);
}

TEST(DistributionsTest, ChebyshevEllipse3Stretched) {
	test_distribution("chebyshev_ellipse_3_stretched", alfi::dist::Type::CHEBYSHEV_ELLIPSE_3_STRETCHED, 1e-13);
}

TEST(DistributionsTest, ChebyshevEllipse4) {
	test_distribution("chebyshev_ellipse_4", alfi::dist::Type::CHEBYSHEV_ELLIPSE_4, 1e-14);
}

TEST(DistributionsTest, ChebyshevEllipse4Stretched) {
	test_distribution("chebyshev_ellipse_4_stretched", alfi::dist::Type::CHEBYSHEV_ELLIPSE_4_STRETCHED, 1e-13);
}

TEST(DistributionsTest, Logistic) {
	test_distribution("logistic", alfi::dist::Type::LOGISTIC, 1e-13);
}

TEST(DistributionsTest, LogisticStretched) {
	test_distribution("logistic_stretched", alfi::dist::Type::LOGISTIC_STRETCHED, 1e-11);
}

TEST(DistributionsTest, Erf) {
	test_distribution("erf", alfi::dist::Type::ERF, 1e-14);
}

TEST(DistributionsTest, ErfStretched) {
	test_distribution("erf_stretched", alfi::dist::Type::ERF_STRETCHED, 1e-11);
}