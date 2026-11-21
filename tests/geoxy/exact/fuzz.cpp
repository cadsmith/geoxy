/***************************************************************************************************
 *  @file      tests/geoxy/exact/fuzz.cpp
 *  @brief     Every input of the fuzz corpus, and every one the fuzzer found, run again.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "../../fuzz/oracle.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <string>
#include <vector>

namespace Geoxy::Fuzz
{
  namespace
  {
    /**
     *  The folder of the fuzz targets, their corpus and what they found.
     */
    auto Inputs() -> std::filesystem::path
    {
      return GEOXY_FUZZ_DIR;
    }

    /**
     *  The files of a folder, by name, none when there is no folder.
     */
    auto Files(const std::filesystem::path &folder) -> std::vector<std::filesystem::path>
    {
      std::vector<std::filesystem::path> files;
      if (!std::filesystem::is_directory(folder))
      {
        return files;
      }
      for (const auto &entry : std::filesystem::directory_iterator{folder})
      {
        if (entry.is_regular_file() && !entry.path().filename().string().starts_with('.'))
        {
          files.push_back(entry.path());
        }
      }
      std::ranges::sort(files);
      return files;
    }

    /**
     *  The bytes of a file.
     */
    auto Bytes(const std::filesystem::path &path) -> std::vector<std::uint8_t>
    {
      std::ifstream     file{path, std::ios::binary};
      const std::string text{std::istreambuf_iterator<char>{file}, {}};
      return {text.begin(), text.end()};
    }

    /**
     *  Runs every input of a folder of each target, checking it breaks nothing.
     */
    void Replay(const std::filesystem::path &folder)
    {
      for (const Named &target : Targets)
      {
        for (const auto &path : Files(folder / target.Name))
        {
          CAPTURE(path.string());
          CHECK(Broken(target.Which, Bytes(path)) == "");
        }
      }
    }
  } // namespace

  TEST_CASE("fuzz: every input of the seed corpus keeps the original's signs and the laws")
  {
    for (const Named &target : Targets)
    {
      CAPTURE(target.Name);
      CHECK_FALSE(Files(Inputs() / "corpus" / target.Name).empty());
    }
    Replay(Inputs() / "corpus");
  }

  TEST_CASE("fuzz: every input the fuzzer found keeps the original's signs and the laws")
  {
    Replay(Inputs() / "found");
  }

  TEST_CASE("fuzz: an empty input draws zeros, and breaks nothing")
  {
    for (const Named &target : Targets)
    {
      CAPTURE(target.Name);
      CHECK(Broken(target.Which, {}) == "");
    }
  }
} // namespace Geoxy::Fuzz
