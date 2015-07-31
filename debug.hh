#ifndef __SIM_DEBUG_H__
#define __SIM_DEBUG_H__


/* Code for printing DEBUG statements.
 * Adapted from http://stackoverflow.com/questions/30073287/c-cout-auto-separator
 */

#include <iostream>
#include <sstream>

struct set_some_separator {
  set_some_separator(const char* sep) : _sep(sep)
  { };

  template <typename T>
    set_some_separator& operator<<(const T& v)
  {
      _str << v << _sep;
        return *this;
    }

  friend
  std::ostream& operator<<(std::ostream& os, const set_some_separator& s)
  { return os << s._str.str(); }

  const char* _sep;
  std::ostringstream _str;
};



#define DEBUG(x) do {							\
    std::cout << "DEBUG " << __FILE__ << " " << __FUNCTION__		\
	      << " " << __LINE__ << " "					\
	      << #x << ": " << (set_some_separator(" ") << x) << std::endl; } \
  while(0)

#define DEBUG_IF(cond, x)			\
  do { if (cond) { DEBUG(x); } } while(0)


#endif
