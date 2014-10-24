#pragma once
#include "typedefs.h"

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

template <typename BitString> inline
bool bitcheck(BitString value, size_t digit)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  return (0x1 & (value >> digit) ) != 0x0;
}

template <typename BitString, typename ...Args> inline
bool bitcheck(BitString value, size_t digit, Args ... args)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  return (value >> digit) & bitcheck(value, args...);
}


template <typename BitString> inline
BitString bitset(BitString value, size_t digit, bool bitvalue)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  BitString mask = (BitString(0x1) << digit);
  value &= ~mask; // set to zero
  mask = bitvalue ? 0x1 : 0x0;
  value |= mask << digit;
  return value;
}

template <> inline
size_t bitcount<uint64_t>(uint64_t val, uint64_t mask)
{
  static const uint64_t S[] = {1, 2, 4, 8, 16, 32}; // Magic Binary Numbers
  static const uint64_t B[] = {0x5555555555555555L, 0x3333333333333333L, 0x0F0F0F0F0F0F0F0FL,
                               0x00FF00FF00FF00FFL, 0x0000FFFF0000FFFFL, 0x00000000FFFFFFFFL};
  uint64_t cnt = 0;
  val = val & mask;
  cnt = val - ((val >> 1) & B[0]);
  cnt = ((cnt >> S[1]) & B[1]) + (cnt & B[1]);
  cnt = ((cnt >> S[2]) + cnt) & B[2];
  cnt = ((cnt >> S[3]) + cnt) & B[3];
  cnt = ((cnt >> S[4]) + cnt) & B[4];
  cnt = ((cnt >> S[5]) + cnt) & B[5];
  return (size_t) cnt;
}

template <> inline
size_t bitcount<uint32_t>(uint32_t val, uint32_t mask)
{
  static const uint32_t S[] = {1, 2, 4, 8, 16}; // Magic Binary Numbers
  static const uint32_t B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFFL};
  uint32_t cnt = 0;
  val = val & mask;
  cnt = val - ((val >> 1) & B[0]);
  cnt = ((cnt >> S[1]) & B[1]) + (cnt & B[1]);
  cnt = ((cnt >> S[2]) + cnt) & B[2];
  cnt = ((cnt >> S[3]) + cnt) & B[3];
  cnt = ((cnt >> S[4]) + cnt) & B[4];
  return (size_t) cnt;
}

template <typename BitString> inline
BitString makemask(size_t first, size_t last)
{
  BitString res = ~BitString(0x0);
  size_t n = sizeof(BitString)*8;
  if (last < first) { std::swap(first, last); }
  
  res = (res << (n-last-1)) >> (n-last-1);
  res = (res >> first) << first;
  return res;
}

template<typename BitString> inline
size_t bitwsum(BitString val, BitString mask)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  size_t sum = 0;
  BitString v = (val & mask)>>1;
  size_t digit = 1;
  while(v) {
    sum += digit * (0x1 & v);
    v >>= 1;
    digit++;
  }
  return sum;
}


template<typename BitString> inline
BitString bitcompress(BitString val, BitString mask)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  BitString result = 0;
  size_t cnt = 0;

  BitString v = val, m = mask;
  while(v && m) {
    if (0x1 & m) {
      result |= BitString(0x1 & v) << cnt;
      cnt++;
    }
    v >>= 1;
    m >>= 1;
  }
  return result;
}

template <typename BitString> inline
BitString bitexpand(BitString val, BitString mask)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  BitString result = 0;
  size_t cnt = 0;
  BitString v = val, m = mask;
  while(m && v) {
    if (0x1 & m) {
      result |= BitString(0x1 & v) << cnt;
      v >>= 1;
    }
    m >>= 1;
    cnt++;
  }
  return result;
}

template <typename BitString> inline
BitString bitswap(BitString val, size_t i1, size_t i2)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  bool b1 = (0x1 & (val >> i1)) != 0;
  bool b2 = (0x1 & (val >> i2)) != 0;
  val = bitset(val, i1, b2);
  val = bitset(val, i2, b1);
  return val;
}



template <typename BitString> inline
size_t partitioncount(BitString val, BitString mask)
{
  static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
  BitString v = val, m = mask;
  size_t ones = 0, count = 0;
  while(v) {
    if(0x1 & v) {
      switch (0x1 & m) {
        case 0x0: count += ones;  break;
        case 0x1: ones++;  break;
      }
    }
    v >>= 1;
    m >>= 1;
  } // i
  return count;
}



namespace naive
{

    template <typename BitString> inline
        BitString makemask(size_t first, size_t last)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            if (last < first) { std::swap(first, last); }
            BitString val = 0x0;
            assert(0 <= first && first <= last && last < sizeof(BitString) * 8);
            BitString one = 0x1;
            for (int i = first; i <= last ; ++i) {
                val |= one << i;
            }
            return val;
        }

    template <typename BitString> inline
        size_t bitcount(BitString val, BitString mask)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            size_t cnt = 0;
            val &= mask;
            for (int i = 0 ; i < sizeof(BitString) * 8 ; ++i) {
                cnt += 0x1 & (val >> i);
            }
            return cnt;
        }

    template<typename BitString> inline
        size_t bitwsum(BitString val, BitString mask)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            size_t sum = 0;
            val = val & mask;
            for (int i = 1 ; i < sizeof(BitString)*8 ; ++i) {
                sum += i * (0x1 & (val>>i) ); 
            }
            return sum;
        }

    template<typename BitString> inline
        BitString bitcompress(BitString val, BitString mask)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            BitString result = 0;
            int cnt = 0;
            for (int i = 0 ; i < sizeof(BitString)*8 ; ++i) {
                if ( 0x1 & (mask>>i) ) {
                    result |= BitString(0x1 & (val>>i)) << cnt;
                    cnt++;
                }
            }
            return result;
        }

    template <typename BitString> inline
        BitString bitexpand(BitString val, BitString mask)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            BitString result = 0x0;
            int cnt = 0;
            for (int i = 0 ; i < sizeof(BitString)*8 ; ++i) {
                if (0x1 & (mask>>i)) {
                    result |= BitString(0x1 & (val>>cnt)) << i;
                    cnt++;
                }
            }
            return result;
        }


    template <typename BitString> inline
        size_t partitioncount(BitString val, BitString mask)
        {
            static_assert(std::is_unsigned<BitString>::value, "BitString should be unsigned");
            size_t count = 0;
            for (int i = 0 ; i < sizeof(BitString)*8 ; ++i) {
                if (0x1 & (~mask >> i) & (val >> i)) {
                    for (int j = 0 ; j < i ; ++j) {
                        if (0x1 & ((mask) >> j) & (val>>j)) {
                            count++;
                        }
                    }
                }
            }
            return count;
        }




} // namespace naive



} // namespace bitbox
} // namespace kore
