#pragma once
#include "../typedefs.h"

//
// Collection of bit operations.
//

namespace kore {
namespace bitbox {

template <typename BitString> inline
bool bitcheck(BitString value, size_t digit);

template <typename BitString> inline
BitString bitset(BitString value, size_t digit, bool bitvalue=true);

template <typename BitString> inline
size_t bitcount(BitString val, BitString mask = ~(BitString(0)));

template <typename BitString> inline
size_t bitwsum(BitString val, BitString mask = ~(BitString(0)));

template <typename BitString> inline
BitString bitcompress(BitString val, BitString mask);

template <typename BitString> inline
BitString bitexpand(BitString val, BitString mask);

template <typename BitString> inline
BitString makemask(size_t first, size_t last);

template <typename BitString> inline
size_t partitioncount(BitString val, BitString mask);

} // namespace bitbox
} // namespace kore

#include "bitbox_impl.h"
