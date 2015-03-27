#pragma once

namespace kore {

#ifdef USE_CXX11
template <typename T, typename U>
using is_same = std::is_same<T, U>;
#else
template <typename A, typename B>
struct is_same { static const bool value = false; };
template <typename A>
struct is_same<A, A> { static const bool value = true; };
#endif

} // namespace kore
