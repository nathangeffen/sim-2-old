#include "test.hh"

using namespace tst;

bool
TestSeries::test(bool expr,
		 const char *description,
		 const char *file,
		 const int line_num)
{
  ++tests_;
  if (expr) {
    ++successes_;
    return true;
  } else {
    ++failures_;
    if (verbose_)
      std::clog << "FAIL:\t" << file << "\t" << line_num << "\t" << description
		<< std::endl;
    return false;
  }
}
