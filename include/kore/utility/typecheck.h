#pragma once

namespace kore {

template <typename A, typename B>
struct is_same { static const bool value = false; };

template <typename A>
struct is_same<A, A> { static const bool value = true; };

template<bool>
struct static_assertion;

template<>
struct static_assertion<true> {};


} // namespace kore
