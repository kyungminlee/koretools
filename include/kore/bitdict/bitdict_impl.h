#pragma once

namespace kore
{
  namespace bitdict
  {
    template <typename BitString, typename Index, typename StringInfo>
    template <typename Reducer> inline
      BitDictionary<BitString, Index, StringInfo>::BitDictionary(Index length, Reducer&& reduce)
      : length_(length), reduce_(reduce)
      , indices_(), words_(), infos_()
    {
      const BitString ONE = 0x1;
      BitString max_word = (ONE << length);
      for (BitString word = 0; word < max_word; ++word) {
        auto reduced = reduce_(word);
        if (!std::get<0>(reduced)) { continue; }

        assert(indices_.find(std::get<1>(reduced)) == indices_.end());

        indices_[std::get<1>(reduced)] = words_.size();
        words_.push_back(std::get<1>(reduced));
        infos_.push_back(std::get<2>(reduced));
      }
    }

    template <typename BitString, typename Index, typename StringInfo> inline
      Index BitDictionary<BitString, Index, StringInfo>::index(BitString wrd) const
    {
      auto reduced = reduce_(wrd);
      if (!std::get<0>(reduced)) { return -1; }
      return indices_.at(std::get<1>(reduced));
    }

    template <typename BitString, typename Index, typename StringInfo> inline
      BitString BitDictionary<BitString, Index, StringInfo>::word(Index idx) const
    {
      assert(0 <= idx && idx < words_.size());
      return words_[idx];
    }

    template <typename BitString, typename Index, typename StringInfo> inline
      const StringInfo& BitDictionary<BitString, Index, StringInfo>::info(Index idx) const
    {
      assert(0 <= idx && idx < words_.size());
      return infos_[idx];
    }

  } // namespace dictionary
} // namespace kore
