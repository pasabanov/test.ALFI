#include <ALFI/ratf.h>
#include <ALFI/dist.h>

#include "../test_utils.h"

TEST(RationalFunctionsTest, Evaluation) {
	alfi::ratf::RationalFunction<> rf({}, {});
	for (const auto x : alfi::dist::uniform(11, -2.0, 2.0)) {
		EXPECT_TRUE(std::isnan(alfi::ratf::val(rf, x)));
	}
	rf = {{1, 0}, {1}};
	for (const auto x : alfi::dist::uniform(11, -2.0, 2.0)) {
		EXPECT_EQ(alfi::ratf::val(rf, x), x);
	}
	rf = {{2, 3, 4}, {1, 2, 2}};
	EXPECT_EQ(alfi::ratf::val(rf, -2.0), 3.0);
	EXPECT_EQ(alfi::ratf::val(rf, -1.0), 3.0);
	EXPECT_EQ(alfi::ratf::val(rf, 0.0), 2.0);
	EXPECT_EQ(alfi::ratf::val(rf, 1.0), 9.0/5.0);
	EXPECT_EQ(alfi::ratf::val(rf, 2.0), 9.0/5.0);
}

TEST(RationalFunctionsTest, Pade) {
	// [0/0] empty array (effectively 0)
	auto rf = alfi::ratf::pade({}, 0, 0);
	expect_eq(rf.first, {0});
	expect_eq(rf.second, {1});
	// [1/1] empty array (effectively 0)
	rf = alfi::ratf::pade({}, 1, 1);
	expect_eq(rf.first, {0});
	expect_eq(rf.second, {1});
	// [1/1] 10
	rf = alfi::ratf::pade({10}, 1, 1);
	expect_eq(rf.first, {10});
	expect_eq(rf.second, {1});
	// [2/2] x^4 + 2x^3 + 3x^2 + 4x + 5
	rf = alfi::ratf::pade({1, 2, 3, 4, 5}, 2, 2);
	expect_eq(rf.first, {-6, 5});
	expect_eq(rf.second, {1, -2, 1});
	// [2/2] x^4 - 2x^3 + 3x^2 - 4x + 5
	rf = alfi::ratf::pade({1, -2, 3, -4, 5}, 2, 2);
	expect_eq(rf.first, {6, 5});
	expect_eq(rf.second, {1, 2, 1});
	// [2/2] exp(x)
	rf = alfi::ratf::pade({1.0/2/3/4, 1.0/2/3, 1.0/2, 1, 1}, 2, 2);
	expect_eq(rf.first, {0.25, 1.5, 3});
	expect_eq(rf.second, {0.25, -1.5, 3});
	// [2/2] sin(x)
	rf = alfi::ratf::pade({1.0/2/3/4, 0, -1.0/2, 0, 1}, 2, 2);
	expect_eq(rf.first, {-5.0/12, 0, 1});
	expect_eq(rf.second, {1.0/12, 0, 1});
	// [2/2] x^5
	rf = alfi::ratf::pade({1, 0, 0, 0, 0, 0}, 2, 2);
	expect_eq(rf.first, {0});
	expect_eq(rf.second, {-1});
	// [2/2] x^4
	rf = alfi::ratf::pade({1, 0, 0, 0, 0}, 2, 2);
	expect_eq(rf.first, {});
	expect_eq(rf.second, {});
}