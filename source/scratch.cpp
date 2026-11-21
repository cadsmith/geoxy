/***************************************************************************************************
 *  @file      source/scratch.cpp
 *  @brief     The arena each thread keeps for what calls work with.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "scratch.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <new>
#include <ranges>
#include <span>

namespace Geoxy
{
  namespace
  {
    constexpr std::align_val_t Alignment = std::align_val_t{16};  /**< Every stretch's alignment. */
    constexpr std::size_t      Initial   = std::size_t{1} << 16U; /**< The first stretch's size.  */
    constexpr std::size_t      Lead      = 16; /**< The bytes before a stretch's first block.     */

    /**
     *  A stretch of memory from the system.
     *
     *  @param  size Its size.
     *  @return      Its first byte.
     */
    auto Asked(std::size_t size) -> std::byte *
    {
      return static_cast<std::byte *>(::operator new(size, Alignment));
    }

    /**
     *  A stretch given back to the system.
     *
     *  @param start Its first byte.
     */
    void Freed(std::byte *start) noexcept
    {
      ::operator delete(start, Alignment);
    }
  } // namespace

  Arena::~Arena()
  {
    for (const std::span<std::byte> stretch : std::span{m_Stretches}.first(m_Count))
    {
      Freed(stretch.data());
    }
  }

  auto Arena::Local() -> Arena &
  {
    thread_local Arena arena;
    return arena;
  }

  auto Arena::Grown(std::size_t size) -> void *
  {
    const std::size_t last  = m_Count == 0 ? 0 : m_Stretches.at(m_Count - 1).size();
    const std::size_t grown = std::max({size + Lead, 2 * last, Initial});
    if (grown > m_Ceiling - std::min(m_Kept, m_Ceiling))
    {
      return Asked(size);
    }
    const std::span whole   = {Asked(grown), grown};
    m_Stretches.at(m_Count) = whole;
    ++m_Count;
    m_Kept += grown;
    m_Room  = whole.subspan(Lead);
    m_Used  = size;
    ++m_Held;
    return m_Room.data();
  }

  auto Arena::Loose(void *room) noexcept -> bool
  {
    const auto start = std::bit_cast<std::uintptr_t>(room);
    const bool held  = std::ranges::any_of(std::span{m_Stretches}.first(m_Count),
                                           [start](std::span<std::byte> stretch) -> bool
                                           {
                                            const auto first =
                                                std::bit_cast<std::uintptr_t>(stretch.data());
                                            return first <= start && start - first < stretch.size();
                                           });
    if (!held)
    {
      Freed(static_cast<std::byte *>(room));
    }
    return !held;
  }

  auto Scratchpad::Shared() -> Scratchpad &
  {
    static Scratchpad shared;
    return shared;
  }

  auto Scratchpad::do_allocate(std::size_t bytes, std::size_t /*alignment*/) -> void *
  {
    return Arena::Local().Take(bytes);
  }

  void Scratchpad::do_deallocate(void *room, std::size_t bytes, std::size_t /*alignment*/)
  {
    Arena::Local().Give(room, bytes);
  }

  auto Scratchpad::do_is_equal(const std::pmr::memory_resource &other) const noexcept -> bool
  {
    return this == &other;
  }

  void Arena::Emptied() noexcept
  {
    const std::span<std::span<std::byte>> stretches = std::span{m_Stretches}.first(m_Count);
    const auto kept = std::ranges::find_if(stretches | std::views::reverse,
                                           [this](std::span<std::byte> stretch) -> bool
                                           { return stretch.size() <= m_Ceiling / 16; });
    const std::span<std::byte> keep =
        kept == (stretches | std::views::reverse).end() ? std::span<std::byte>{} : *kept;
    for (const std::span<std::byte> stretch : stretches)
    {
      if (stretch.data() != keep.data())
      {
        Freed(stretch.data());
      }
    }
    m_Count        = keep.empty() ? 0 : 1;
    m_Stretches[0] = keep;
    m_Kept         = keep.size();
    m_Room         = keep.empty() ? keep : keep.subspan(Lead);
    m_Used         = 0;
  }
} // namespace Geoxy
