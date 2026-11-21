/***************************************************************************************************
 *  @file      tests/suites/cases.hpp
 *  @brief     The cases of the suites, as their files hold them: two curves, a tolerance, the hits
 *             and shared stretches expected.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_SUITES_CASES_HPP
#define GEOXY_TESTS_SUITES_CASES_HPP

#include "support/text.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <vector>

namespace Geoxy::Suites
{
  /**
   *  A hit a case expects: where, how the curves meet there, if said, and how far off geoxy's may
   *  be, which along a tangency is the band where the curves stay within the tolerance.
   */
  struct Expected
  {
    Point       At    = {};  /**< Where.                                          */
    std::string Kind  = {};  /**< Cross, Touch or End; empty if not said.         */
    double      Slack = 0.0; /**< How far off a hit may be; 0 for the case's own. */
  };

  /**
   *  A case: two curves, the tolerance to meet them by, and what they should give.
   */
  struct Case
  {
    std::string                       Name     = {};    /**< Where it is from, and its title. */
    double                            Within   = 0.0;   /**< The tolerance.                   */
    bool                              Relative = false; /**< Whether the tolerance scales.    */
    bool                              Disabled = false; /**< Whether its source skips it.     */
    std::string                       Limit    = {};    /**< Why geoxy cannot hold it, if so. */
    std::optional<Curve>              First    = {};    /**< The first curve.                 */
    std::optional<Curve>              Second   = {};    /**< The second.                      */
    std::vector<Expected>             Points   = {};    /**< The hits expected.               */
    std::vector<std::array<Point, 2>> Shared   = {};    /**< The ends of each shared arc.     */
  };

  /**
   *  The cases in a file.
   *
   *  @param  file The file.
   *  @return      Its cases.
   */
  inline auto Read(const std::filesystem::path &file) -> std::vector<Case>
  {
    std::vector<Case> cases;
    std::ifstream     input{file};
    Case              current;
    for (std::string line; std::getline(input, line);)
    {
      std::istringstream words{line};
      std::string        word;
      words >> word;
      const std::string rest = line.substr(std::min(line.size(), word.size() + 1));
      if (word == "case")
      {
        current      = {};
        current.Name = rest;
      }
      else if (word == "tolerance")
      {
        words >> current.Within;
        current.Relative = rest.contains("relative");
      }
      else if (word == "disabled")
      {
        current.Disabled = true;
      }
      else if (word == "limit")
      {
        current.Limit = rest;
      }
      else if (word == "first")
      {
        current.First = FromText(rest);
      }
      else if (word == "second")
      {
        current.Second = FromText(rest);
      }
      else if (word == "point")
      {
        Expected point;
        words >> point.At.X >> point.At.Y;
        current.Points.push_back(point);
      }
      else if (word == "hit")
      {
        Expected hit;
        words >> hit.At.X >> hit.At.Y >> hit.Kind >> hit.Slack;
        current.Points.push_back(hit);
      }
      else if (word == "shared")
      {
        std::array<Point, 2> ends;
        for (double *number : {&ends[0].X, &ends[0].Y, &ends[1].X, &ends[1].Y})
        {
          std::string text;
          words >> text;
          *number = std::stod(text);
        }
        current.Shared.push_back(ends);
      }
      else if (word == "end")
      {
        cases.push_back(current);
      }
    }
    return cases;
  }

  /**
   *  Every case file under some roots, in a folder named cases; a root that is not there gives
   *  none.
   *
   *  @param  roots The roots.
   *  @return       The files.
   */
  inline auto Files(std::span<const std::filesystem::path> roots)
      -> std::vector<std::filesystem::path>
  {
    std::vector<std::filesystem::path> files;
    for (const std::filesystem::path &root : roots)
    {
      if (!std::filesystem::is_directory(root))
      {
        continue;
      }
      for (const auto &entry : std::filesystem::recursive_directory_iterator{root})
      {
        if (entry.path().extension() == ".txt" && entry.path().parent_path().filename() == "cases")
        {
          files.push_back(entry.path());
        }
      }
    }
    std::ranges::sort(files);
    return files;
  }
} // namespace Geoxy::Suites

#endif /* GEOXY_TESTS_SUITES_CASES_HPP */
