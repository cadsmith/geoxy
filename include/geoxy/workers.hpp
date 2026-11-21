/***************************************************************************************************
 *  @file      geoxy/workers.hpp
 *  @brief     Threads kept by a caller to work a call on at once, the calling thread among them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WORKERS_HPP
#define GEOXY_WORKERS_HPP

#include "geoxy/api.hpp"

#include <cstddef>
#include <memory>

namespace Geoxy
{
  /**
   *  Threads a caller keeps, to work a call that takes them on all at once: the calling thread, and
   *  as many more as asked, which wait between calls and are woken for each.
   *
   *  A pool is the caller's own, made once and handed to every call that is to use it, so that no
   *  call pays for threads to start: there is none of geoxy's own anywhere. Each thread works in
   *  its own memory, kept between calls. What a call gives is the same, to the bit and in its
   *  order, whatever the pool, as working alone gives it.
   *
   *  A pool works one call at a time. A call that finds it working another, one made at once on
   *  another thread or from within what the call tells, works alone on the calling thread, which
   *  gives the same, and waits for no one. An exception thrown by what a call tells is thrown on to
   *  the caller once the pool's threads are done with the call; one on any of those threads ends
   *  the program, as the standard's parallel algorithms do.
   */
  class Workers
  {
  public:
    /**
     *  A pool of threads, the calling thread counted among them.
     *
     *  @param threads How many, the calling thread among them: as many as the hardware runs at once
     *                 for none, one for the calling thread alone.
     */
    GEOXY_API explicit Workers(std::size_t threads = 0);

    /**
     *  No copy: a thread is one pool's.
     */
    Workers(const Workers &) = delete;

    /**
     *  No move: the threads hold where their pool is.
     */
    Workers(Workers &&) = delete;

    /**
     *  No copy: a thread is one pool's.
     *
     *  @return Nothing.
     */
    auto operator=(const Workers &) -> Workers & = delete;

    /**
     *  No move: the threads hold where their pool is.
     *
     *  @return Nothing.
     */
    auto operator=(Workers &&) -> Workers & = delete;

    /**
     *  The threads stopped and let go, once each has ended what it was working.
     */
    GEOXY_API ~Workers();

    /**
     *  How many threads it works a call on, the calling thread among them.
     *
     *  @return The count, one at the least.
     */
    [[nodiscard]] GEOXY_API auto Size() const -> std::size_t;

  private:
    friend struct Pooled;

    struct Shared;

    std::unique_ptr<Shared> m_Shared; /**< The threads, and what they share. */
  };
} // namespace Geoxy

#endif /* GEOXY_WORKERS_HPP */
