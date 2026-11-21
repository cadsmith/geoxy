/***************************************************************************************************
 *  @file      source/unions.hpp
 *  @brief     Sets of indices joined one pair at a time, each stood for by its least member.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_UNIONS_HPP
#define GEOXY_UNIONS_HPP

#include "scratch.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>

namespace Geoxy
{
  /**
   *  Sets of indices joined one pair at a time, each stood for by its least member, so that what
   *  stands for a set does not hang on the order the pairs were joined in. The parents lie in the
   *  thread's scratch memory: the sets are made, used and dropped on one thread.
   */
  class Unions
  {
  public:
    /**
     *  Each of some indices a set of its own.
     *
     *  @param count How many.
     */
    explicit Unions(std::size_t count) : m_Parents(count)
    {
      std::ranges::iota(m_Parents, std::size_t{0});
    }

    /**
     *  The index that stands for an index's set.
     *
     *  @param  member The index.
     *  @return        The least index of its set.
     */
    auto Root(std::size_t member) -> std::size_t
    {
      std::size_t root = member;
      while (m_Parents[root] != root)
      {
        root = m_Parents[root];
      }
      while (m_Parents[member] != root)
      {
        member = std::exchange(m_Parents[member], root);
      }
      return root;
    }

    /**
     *  Two indices' sets made one.
     *
     *  @param first  One index.
     *  @param second The other.
     */
    void Join(std::size_t first, std::size_t second)
    {
      const std::size_t one           = Root(first);
      const std::size_t other         = Root(second);
      m_Parents[std::max(one, other)] = std::min(one, other);
    }

  private:
    List<std::size_t> m_Parents = {}; /**< Each index's parent. */
  };
} // namespace Geoxy

#endif /* GEOXY_UNIONS_HPP */
