/***************************************************************************************************
 *  @file      source/workers.cpp
 *  @brief     A pool's threads, woken for each job by a round that moves on, the round ended by the
 *             last of them done with it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/workers.hpp"

#include "pooled.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <stop_token>
#include <thread>
#include <vector>

namespace Geoxy
{
  /**
   *  What a pool's threads share: the round they are on, the job of that round, and how many of
   *  them are done with it.
   *
   *  The caller sets the job, then moves the round on; a thread that sees the round move works the
   *  job, then counts itself done, as the caller does; the last to, of all of them, ends the round,
   *  which the caller waits for. So no thread starts a round before every one has ended the last,
   *  and the job is not set again while one reads it. Every count is an exchange that acquires what
   *  was counted before it and releases what its thread did, so that the round's end follows all
   *  each thread did in it.
   */
  struct Workers::Shared
  {
    /**
     *  What a pool of some threads shares, with no thread started yet.
     *
     *  @param size How many threads, the calling thread among them.
     */
    explicit Shared(std::size_t size) : Size(size) {}

    Shared(const Shared &)                     = delete;
    Shared(Shared &&)                          = delete;
    auto operator=(const Shared &) -> Shared & = delete;
    auto operator=(Shared &&) -> Shared &      = delete;

    /**
     *  Every thread told to stop and woken, then let go once it has.
     */
    ~Shared()
    {
      for (std::jthread &thread : Threads)
      {
        thread.request_stop();
      }
      Round.fetch_add(1, std::memory_order_release);
      Round.notify_all();
      Threads.clear();
    }

    std::size_t                Size  = 1;        /**< The threads, the caller among them. */
    std::atomic<std::uint64_t> Round = 0;        /**< How many rounds have started.       */
    std::atomic<std::uint64_t> Done  = 0;        /**< How many threads were done with one,
                                                      over every round.                     */
    std::atomic<std::uint64_t> Ended = 0;        /**< How many rounds every thread is done
                                                      with.                                 */
    std::atomic<bool>         Busy    = false;   /**< Whether a call holds the pool. */
    const Pooled::Job        *Job     = nullptr; /**< The job of the round.          */
    std::vector<std::jthread> Threads = {};      /**< The threads but the caller, joined last
                                                      as they are let go first.             */
  };

  namespace
  {
    /**
     *  The round a pool is on, once it has moved past one.
     *
     *  @param  round The pool's round.
     *  @param  seen  The round seen last.
     *  @return       The round now.
     */
    auto Next(const std::atomic<std::uint64_t> &round, std::uint64_t seen) -> std::uint64_t
    {
      round.wait(seen, std::memory_order_acquire);
      return round.load(std::memory_order_acquire);
    }
  } // namespace

  Workers::Workers(std::size_t threads)
      : m_Shared(std::make_unique<Shared>(
            threads == 0 ? std::max<std::size_t>(1, std::thread::hardware_concurrency()) : threads))
  {
    m_Shared->Threads.reserve(m_Shared->Size - 1);
    for (std::size_t slot = 1; slot < m_Shared->Size; ++slot)
    {
      m_Shared->Threads.emplace_back(&Pooled::Worked, std::ref(*m_Shared), slot);
    }
  }

  Workers::~Workers() = default;

  auto Workers::Size() const -> std::size_t
  {
    return m_Shared->Size;
  }

  void Pooled::Worked(const std::stop_token &stop, Workers::Shared &shared, std::size_t slot)
  {
    for (std::uint64_t seen = Next(shared.Round, 0); !stop.stop_requested();
         seen               = Next(shared.Round, seen))
    {
      (*shared.Job)(slot);
      Ended(shared, seen);
    }
  }

  void Pooled::Ended(Workers::Shared &shared, std::uint64_t round)
  {
    if (shared.Done.fetch_add(1, std::memory_order_acq_rel) + 1 == round * shared.Size)
    {
      shared.Ended.store(round, std::memory_order_release);
      shared.Ended.notify_all();
    }
  }

  auto Pooled::Run(Workers &pool, const Job &job) -> bool
  {
    Workers::Shared &shared = *pool.m_Shared;
    if (shared.Busy.exchange(true, std::memory_order_acquire))
    {
      return false;
    }
    shared.Job                = &job;
    const std::uint64_t round = shared.Round.fetch_add(1, std::memory_order_release) + 1;
    shared.Round.notify_all();
    std::exception_ptr thrown;
    try
    {
      job(0);
    }
    catch (...)
    {
      thrown = std::current_exception();
    }
    Ended(shared, round);
    shared.Ended.wait(round - 1, std::memory_order_acquire);
    shared.Busy.store(false, std::memory_order_release);
    if (thrown)
    {
      std::rethrow_exception(thrown);
    }
    return true;
  }
} // namespace Geoxy
