#ifndef HANDY_CHRONO_H
#define HANDY_CHRONO_H

#include <chrono>

using chrono_clock = std::chrono::time_point<std::chrono::high_resolution_clock>;

inline chrono_clock now()
{
  return std::chrono::high_resolution_clock::now();
}

template <typename T>
inline auto milliseconds(T t)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(t).count();
}

template <typename T>
inline auto microseconds(T t)
{
  return std::chrono::duration_cast<std::chrono::microseconds>(t).count();
}

template <typename T>
inline auto nanoseconds(T t)
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(t).count();
}

inline double duration(chrono_clock start, chrono_clock end)
{
  return milliseconds(end - start) / 1e3;
}

inline double higher_precision_duration(chrono_clock start, chrono_clock end)
{
  return nanoseconds(end - start) / 1e9;
}

#endif // HANDY_CHRONO_H
