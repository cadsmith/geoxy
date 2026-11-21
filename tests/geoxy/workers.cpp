/***************************************************************************************************
 *  @file      tests/geoxy/workers.cpp
 *  @brief     Tests that a pool works a job on each of its threads once a call, the calling thread
 *             among them, however many calls it is given, and throws on what the caller's throws.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/workers.hpp"

#include "pooled.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <stdexcept>
#include <thread>

namespace Geoxy
{
  TEST_CASE("workers: as many threads as asked, or as the hardware runs at once")
  {
    const Workers most;
    CHECK(most.Size()       == std::max<std::size_t>(1, std::thread::hardware_concurrency()));
    CHECK(Workers{1}.Size() == 1);
    CHECK(Workers{7}.Size() == 7);
  }

  TEST_CASE("workers: each thread works each job once, the calling thread as the first")
  {
    Workers                         pool{5};
    std::array<std::atomic<int>, 5> worked = {};
    std::atomic<std::thread::id>    first;
    const Pooled::Job               job = [&worked, &first](std::size_t slot) -> void
    {
      worked.at(slot).fetch_add(1);
      if (slot == 0)
      {
        first.store(std::this_thread::get_id());
      }
    };
    for (int round = 0; round < 40; ++round)
    {
      CHECK(Pooled::Run(pool, job));
    }
    CHECK(std::ranges::all_of(worked, [](const std::atomic<int> &count) -> bool
                              { return count.load() == 40; }));
    CHECK(first.load() == std::this_thread::get_id());

    Workers alone{1};
    CHECK(Pooled::Run(alone, job));
    CHECK(worked.at(0).load() == 41);
  }

  TEST_CASE(
      "workers: what the job throws on the calling thread is thrown on once the others are done")
  {
    Workers           pool{4};
    std::atomic<int>  done = 0;
    const Pooled::Job job  = [&done](std::size_t slot) -> void
    {
      if (slot == 0)
      {
        throw std::runtime_error{"thrown"};
      }
      done.fetch_add(1);
    };
    CHECK_THROWS_AS(Pooled::Run(pool, job), std::runtime_error);
    CHECK(done.load() == 3);
    CHECK(Pooled::Run(pool, [&done](std::size_t) -> void { done.fetch_add(1); }));
    CHECK(done.load() == 7);
  }

  TEST_CASE("workers: a pool working one job is not given another, and a pool never given one ends")
  {
    Workers           pool{3};
    bool              inner = true;
    const Pooled::Job job   = [&pool, &inner](std::size_t slot) -> void
    {
      if (slot == 0)
      {
        inner = Pooled::Run(pool, [](std::size_t) -> void {});
      }
    };
    CHECK(Pooled::Run(pool, job));
    CHECK(!inner);
    const Workers idle{6};
    CHECK(idle.Size() == 6);
  }
} // namespace Geoxy
