#pragma once

#include <exception>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <functional>
#include <type_traits>
#include "../typedefs.h"

namespace kore
{
  namespace bitdict
  {
    template <typename _BitString = uint64_t,
      typename _Index = int64_t,
      typename _StringInfo = void*>
    class BitDictionary
    {
    public:
      typedef _BitString BitString;
      typedef _Index Index;
      typedef _StringInfo StringInfo;

      static_assert(sizeof(Index) >= sizeof(BitString),
                    "BitDictionary: size of Index and size of BitString should match");
      static_assert(std::is_unsigned<BitString>::value && std::is_integral<BitString>::value,
                    "BitDictionary: BitString should be unsigned integer");
      static_assert(std::is_signed<Index>::value && std::is_integral<Index>::value,
                    "BitDictionary: Index should be signed integer");

      template <typename Reducer>
      BitDictionary(Index length, Reducer&& reduce);

      Index size() const { return words_.size();}
      Index index(BitString word) const;
      BitString word(Index idx) const;
      const StringInfo& info(Index idx) const;
      std::function<std::tuple<bool, BitString, StringInfo>(BitString)>& reduce() { return reduce_; }
    private:
      Index length_;
      std::function<std::tuple<bool, BitString, StringInfo>(BitString) > reduce_;
      std::unordered_map<BitString, Index> indices_;
      std::vector<BitString> words_;
      std::vector<StringInfo> infos_;
    };

  } // namespace dictionary
} // namespace kore

#include "bitdict_impl.h"
