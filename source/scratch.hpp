/***************************************************************************************************
 *  @file      source/scratch.hpp
 *  @brief     Lists that draw on an arena each thread keeps: room is taken by moving a mark, the
 *             last room taken given back by moving it back, and the whole arena is empty again once
 *             every list drawing on it is gone, so that working asks the system for no memory once
 *             the arena has grown to what a call needs. Past a ceiling, room comes from the system
 *             and goes straight back to it, so no call holds more than the ceiling in the arena.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SCRATCH_HPP
#define GEOXY_SCRATCH_HPP

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <span>
#include <type_traits>
#include <vector>

namespace Geoxy
{
  /**
   *  The memory a thread works in: blocks of it handed out in turn from a mark that moves on, the
   *  latest one handed out taken back by moving the mark back, any other left where it is until
   *  every block is back, when the arena starts over from its first byte. When the room runs out a
   *  stretch twice as large is added, while the stretches together stay within the ceiling; past
   *  it, a block is the system's own, given back the moment it is. Blocks given back out of turn
   *  stay taken until the arena starts over, so without the ceiling a long call would hold ever
   *  more. When the arena starts over it keeps its largest stretch within a sixteenth of the
   *  ceiling, so that a call works in one stretch once the arena has seen its like, and no thread
   *  holds more than that between calls. Each stretch's blocks start 16 bytes into it, so that no
   *  block of another stretch, nor one of the system's, ends where its mark stands.
   */
  class Arena
  {
  public:
    /**
     *  The bytes an arena's stretches hold together, at the most, unless it is given another.
     */
    static constexpr std::size_t Ceiling = std::size_t{1} << 26U;

    /**
     *  An empty arena.
     *
     *  @param ceiling The bytes its stretches may hold together, at the most.
     */
    explicit Arena(std::size_t ceiling = Ceiling) noexcept : m_Ceiling(ceiling) {}

    Arena(const Arena &)                     = delete;
    Arena(Arena &&)                          = delete;
    auto operator=(const Arena &) -> Arena & = delete;
    auto operator=(Arena &&) -> Arena &      = delete;
    ~Arena();

    /**
     *  The arena of the calling thread.
     *
     *  @return The arena.
     */
    static auto Local() -> Arena &;

    /**
     *  Room for some bytes, aligned to 16.
     *
     *  @param  bytes How many.
     *  @return       The room.
     */
    [[nodiscard]] auto Take(std::size_t bytes) -> void *
    {
      const std::size_t size = Rounded(bytes);
      if (size <= m_Room.size() - m_Used)
      {
        void *room  = m_Room.subspan(m_Used, size).data();
        m_Used     += size;
        ++m_Held;
        return room;
      }
      return Grown(size);
    }

    /**
     *  Room given back: the mark moved back where it was the latest taken, the arena started over
     *  where it was the last held.
     *
     *  @param room  The room.
     *  @param bytes How many bytes it was taken for.
     */
    void Give(void *room, std::size_t bytes) noexcept
    {
      const std::size_t size  = Rounded(bytes);
      const auto        start = std::bit_cast<std::uintptr_t>(room);
      const auto        mark  = std::bit_cast<std::uintptr_t>(m_Room.data()) + m_Used;
      if (start + size == mark)
      {
        m_Used -= size;
      }
      else if (Loose(room))
      {
        return;
      }
      --m_Held;
      if (m_Held == 0)
      {
        Emptied();
      }
    }

    /**
     *  The bytes its stretches hold together.
     *
     *  @return How many.
     */
    [[nodiscard]] auto Kept() const noexcept -> std::size_t { return m_Kept; }

  private:
    /**
     *  A byte count rounded up to a whole number of 16-byte steps, one at the least.
     *
     *  @param  bytes The count.
     *  @return       It rounded.
     */
    static constexpr auto Rounded(std::size_t bytes) -> std::size_t
    {
      constexpr std::size_t Step = 16;
      return bytes == 0 ? Step : (bytes + Step - 1) & ~(Step - 1);
    }

    /**
     *  Room taken from a new, larger stretch, where the current one has too little left; or the
     *  system's own, where a new stretch would take the arena past its ceiling.
     *
     *  @param  size How many bytes, rounded.
     *  @return      The room.
     */
    auto Grown(std::size_t size) -> void *;

    /**
     *  Whether room given back out of turn is the system's own, lying in no stretch, and if so
     *  given back to the system.
     *
     *  @param  room The room.
     *  @return      Whether it was.
     */
    auto Loose(void *room) noexcept -> bool;

    /**
     *  The arena started over, every block back: only its largest stretch within a sixteenth of its
     *  ceiling kept, if any is.
     */
    void Emptied() noexcept;

    /**
     *  How many stretches the arena may have: each twice the last, from 64 KiB, so that 48 hold
     *  more than any memory, and any ceiling is reached first.
     */
    static constexpr std::size_t Most = 48;

    std::array<std::span<std::byte>, Most> m_Stretches = {}; /**< Each stretch, whole.            */
    std::size_t                            m_Count     = 0;  /**< How many there are.             */
    std::span<std::byte>                   m_Room      = {}; /**< The current one, past its lead. */
    std::size_t                            m_Used      = 0;  /**< How much of it is taken.        */
    std::size_t                            m_Held      = 0;  /**< How many blocks are out.        */
    std::size_t                            m_Kept      = 0;  /**< The stretches' bytes together.  */
    std::size_t                            m_Ceiling   = 0;  /**< The most they may hold.         */
  };

  /**
   *  An allocator that draws on the arena of the thread that made it, read once when it is made, so
   *  that taking and giving back room never looks for it again. Copies draw on the same arena; a
   *  list is used and dropped on the thread that made it.
   *
   *  @tparam ValueT What it allocates.
   */
  template <typename ValueT>
  class Scratch
  {
  public:
    using value_type      = ValueT;         /**< What it allocates.                     */
    using is_always_equal = std::true_type; /**< Every one draws on its thread's arena. */

    Scratch() = default;

    /**
     *  The allocator of another kind, the same arena.
     *
     *  @param other The allocator.
     */
    template <typename OtherT>
    constexpr explicit(false) Scratch(const Scratch<OtherT> &other) noexcept : m_Arena(other.Seen())
    {
    }

    /**
     *  Room for some values, from the arena.
     *
     *  @param  count How many.
     *  @return       The room.
     */
    [[nodiscard]] auto allocate(std::size_t count) -> ValueT *
    {
      return static_cast<ValueT *>(Drawn().Take(count * sizeof(ValueT)));
    }

    /**
     *  Room given back to the arena.
     *
     *  @param values The room.
     *  @param count  How many values it held.
     */
    void deallocate(ValueT *values, std::size_t count) noexcept
    {
      Drawn().Give(values, count * sizeof(ValueT));
    }

    /**
     *  The arena it has drawn on so far, if any.
     *
     *  @return The arena, or nothing.
     */
    [[nodiscard]] auto Seen() const noexcept -> Arena * { return m_Arena; }

    /**
     *  Whether two allocators are one: always, as each draws on its thread's arena.
     *
     *  @return True.
     */
    [[nodiscard]] friend constexpr auto operator==(const Scratch & /*lhs*/, const Scratch & /*rhs*/)
        -> bool
    {
      return true;
    }

  private:
    /**
     *  The arena it draws on: the calling thread's, looked for once, the first time it is needed.
     *
     *  @return The arena.
     */
    auto Drawn() noexcept -> Arena &
    {
      if (m_Arena == nullptr)
      {
        m_Arena = &Arena::Local();
      }
      return *m_Arena;
    }

    Arena *m_Arena = nullptr; /**< The arena, once looked for. */
  };

  /**
   *  The calling thread's arena as a memory resource, for lists of the standard's kind made, used
   *  and dropped on one thread, as Scratch is for List: room taken from the arena and given back to
   *  it, aligned to 16, as much as anything a list holds asks for.
   */
  class Scratchpad final : public std::pmr::memory_resource
  {
  public:
    /**
     *  The one resource all threads share, each drawing on its own arena.
     *
     *  @return The resource.
     */
    static auto Shared() -> Scratchpad &;

  private:
    /**
     *  Room for some bytes, from the calling thread's arena.
     *
     *  @param  bytes How many.
     *  @return       The room, aligned to 16.
     */
    auto do_allocate(std::size_t bytes, std::size_t /*alignment*/) -> void * override;

    /**
     *  Room given back to the calling thread's arena.
     *
     *  @param room  The room.
     *  @param bytes How many bytes it was taken for.
     */
    void do_deallocate(void *room, std::size_t bytes, std::size_t /*alignment*/) override;

    /**
     *  Whether another resource is this one: room from one can be given back only to it.
     *
     *  @param  other The other.
     *  @return       Whether it is.
     */
    [[nodiscard]] auto do_is_equal(const std::pmr::memory_resource &other) const noexcept
        -> bool override;
  };

  /**
   *  A list drawing on the calling thread's arena.
   *
   *  @tparam ValueT What it holds.
   */
  template <typename ValueT>
  using List = std::vector<ValueT, Scratch<ValueT>>;
} // namespace Geoxy

#endif /* GEOXY_SCRATCH_HPP */
