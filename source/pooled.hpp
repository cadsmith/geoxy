/***************************************************************************************************
 *  @file      source/pooled.hpp
 *  @brief     What the library asks of a pool: to be held for a call, each of its threads working
 *             one job, the calling thread among them, and the call ending once every one is done.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_POOLED_HPP
#define GEOXY_POOLED_HPP

#include "geoxy/workers.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stop_token>

namespace Geoxy
{
  /**
   *  What the library asks of a pool.
   */
  struct Pooled
  {
    /**
     *  A job each thread of a pool works, told which thread it is: the calling thread is 0, the
     *  others 1 and on.
     */
    using Job = std::function<void(std::size_t slot)>;

    /**
     *  A job worked by every thread of a pool at once, the calling thread among them, where the
     *  pool is free: the call returns once every thread is done with it. What the job throws on the
     *  calling thread is thrown on from here once the others are done.
     *
     *  @param  pool The pool.
     *  @param  job  The job.
     *  @return      Whether the pool was free, and the job worked; nothing is done where it was
     *               working another call.
     */
    static auto Run(Workers &pool, const Job &job) -> bool;

    /**
     *  What a thread of a pool does: wait for a job, work it, and wait for the next, until the pool
     *  is let go.
     *
     *  @param stop   Whether the pool is let go.
     *  @param shared What the pool's threads share.
     *  @param slot   Which thread it is, 1 and on.
     */
    static void Worked(const std::stop_token &stop, Workers::Shared &shared, std::size_t slot);

    /**
     *  A thread of a pool counted done with a round, and the round ended where it is the last.
     *
     *  @param shared What the pool's threads share.
     *  @param round  The round.
     */
    static void Ended(Workers::Shared &shared, std::uint64_t round);

    /**
     *  A flag set and told to every thread waiting on it when what it marks is done, however it
     *  ends, so that no thread waits on work that threw.
     */
    class Done
    {
    public:
      /**
       *  A flag to set once a scope ends.
       *
       *  @param flag The flag.
       */
      explicit Done(std::atomic<bool> &flag) noexcept : m_Flag(&flag) {}

      Done(const Done &)                     = delete;
      Done(Done &&)                          = delete;
      auto operator=(const Done &) -> Done & = delete;
      auto operator=(Done &&) -> Done &      = delete;

      /**
       *  The flag set, and every thread waiting on it woken.
       */
      ~Done()
      {
        m_Flag->store(true, std::memory_order_release);
        m_Flag->notify_all();
      }

    private:
      std::atomic<bool> *m_Flag = nullptr; /**< The flag. */
    };
  };
} // namespace Geoxy

#endif /* GEOXY_POOLED_HPP */
