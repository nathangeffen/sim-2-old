#ifndef __TEST_H__
#define __TEST_H__

/**
 * Simple test framework.
 *
 */
#include <iomanip>
#include <iostream>
#include <string>

namespace tst {

  class TestSeries {
  private:
    unsigned tests_ = 0;
    unsigned successes_ = 0;
    unsigned failures_ = 0;
    std::string description_ = 0;
    bool verbose_;
  public:
    TestSeries(const char* description = "",
	       const bool verbose = true) : description_(description),
					    verbose_(verbose) {}
    inline unsigned successes() const { return successes_; }
    inline unsigned failures() const { return failures_; }
    inline unsigned tests() const { return tests_; }
    void set_verbose(const bool verbose) { verbose_ = verbose; }
    bool get_verbose() const { return verbose_; }
    bool test(bool expr,
	      const char *description,
	      const char *file,
	      const int line_num);
    unsigned summary() {
      std::clog << description_ << " tests:\t" << tests_
		<< "\tSuccesses:\t" << successes_
		<< "\tFailures:\t" << failures_ << std::endl;
      return failures_;
    }
  };
}

#define TEST(test_series, expr, desc) \
  test_series.test(expr, desc, __FILE__, __LINE__)

#define TESTCMP(test_series, ex1, cmp, ex2, desc) \
  if ( TEST(test_series, ex1 cmp ex2, desc) == false) {		\
    std::clog << # ex1 << " " << #cmp << " " << # ex2 << std::endl;	\
      std::clog << std::fixed << (ex1) << " " << #cmp << " " \
		<< std::fixed << (ex2) << std::endl;	\
  }

#define TESTEQ(test_series, ex1, ex2, desc) \
  TESTCMP(test_series, ex1, ==, ex2, desc)

#define TESTLT(test_series, ex1, ex2, desc) \
  TESTCMP(test_series, ex1, <, ex2, desc)

#define TESTDBL(test_series, ex1, ex2, desc) \
  TESTCMP(test_series, (ex1) - (ex2), < , 0.00000000000000000001, desc)

#endif
