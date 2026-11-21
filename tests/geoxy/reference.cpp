/***************************************************************************************************
 *  @file      tests/geoxy/reference.cpp
 *  @brief     Offsets held to the answers the reference program gives: shapes drawn at random, each
 *             answer worked out to 120 digits and rounded once, and geoxy's the same to the bit.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/extents.hpp"
#include "geoxy/length.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A case of tests/geoxy/cases/offsets.txt: a shape, what is asked of it, and the answer.
     */
    struct Reference
    {
      std::string           Name     = {};    /**< Its name.                               */
      Polyline              Shape    = {};    /**< The shape.                              */
      bool                  Buffered = false; /**< Whether a region is asked, or one side. */
      Offsetting            Offset   = {};    /**< One side: the distance and the gaps.    */
      Buffering             Buffer   = {};    /**< A region: the distance and how.         */
      Tolerance             Within   = {};    /**< The tolerance it was answered within.   */
      std::vector<Polyline> Answer   = {};    /**< The answer, as geoxy writes it.         */
    };

    /**
     *  A polyline read from a line's words: whether it is closed, how many vertices, and each
     *  vertex's x, y and bulge in hex.
     *
     *  @param  words The words, after the line's first.
     *  @return       The polyline.
     */
    auto Drawn(std::istringstream &words) -> Polyline
    {
      Polyline    polyline;
      int         closed = 0;
      std::size_t count  = 0;
      words >> closed >> count;
      polyline.Closed = closed != 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        std::string x;
        std::string y;
        std::string bulge;
        words >> x >> y >> bulge;
        polyline.Vertices.push_back({.Position = {.X = std::strtod(x.c_str(), nullptr),
                                                  .Y = std::strtod(y.c_str(), nullptr)},
                                     .Bulge    = std::strtod(bulge.c_str(), nullptr)});
      }
      return polyline;
    }

    /**
     *  How a corner's gap is closed, from the suite's word for it.
     *
     *  @param  word Round, Mitre or Chamfer.
     *  @return      The gap kind.
     */
    auto Gapped(const std::string &word) -> Gap
    {
      if (word == "Round")
      {
        return Gap::Fillet;
      }
      return word == "Mitre" ? Gap::Extend : Gap::Chamfer;
    }

    /**
     *  Which way a region is offset, from the suite's word for it.
     *
     *  @param  word Inward, Outward or Both.
     *  @return      The side.
     */
    auto Sided(const std::string &word) -> Side
    {
      if (word == "Inward")
      {
        return Side::Inward;
      }
      return word == "Outward" ? Side::Outward : Side::Both;
    }

    /**
     *  How an open shape's ends are closed, from the suite's word for it.
     *
     *  @param  word Round, Butt or Square.
     *  @return      The cap.
     */
    auto Capped(const std::string &word) -> Cap
    {
      if (word == "Round")
      {
        return Cap::Round;
      }
      return word == "Butt" ? Cap::Flat : Cap::Square;
    }

    /**
     *  What a region's line asks: the distance, the side, the filling, the corners and the ends.
     *
     *  @param  words The words, after the line's first.
     *  @return       The buffering.
     */
    auto Asked(std::istringstream &words) -> Buffering
    {
      std::string distance;
      std::string side;
      std::string fill;
      std::string corner;
      std::string cap;
      words >> distance >> side >> fill >> corner >> cap;
      return {.Distance = std::strtod(distance.c_str(), nullptr),
              .Sides    = Sided(side),
              .Fill     = fill == "EvenOdd" ? Fills::EvenOdd : Fills::NonZero,
              .Gaps     = Gapped(corner),
              .Ends     = Capped(cap)};
    }

    /**
     *  Every case of the file.
     *
     *  @return The cases, in the file's order.
     */
    auto References() -> std::vector<Reference>
    {
      std::vector<Reference> cases;
      std::ifstream          file{
          std::filesystem::path{GEOXY_TESTS_DIR} / "geoxy" / "cases" / "offsets.txt"};
      Reference one;
      for (std::string line; std::getline(file, line);)
      {
        std::istringstream words{line};
        std::string        first;
        words >> first;
        if (first == "case")
        {
          one = Reference{};
          words >> one.Name;
        }
        else if (first == "shape")
        {
          one.Shape = Drawn(words);
        }
        else if (first == "buffer")
        {
          one.Buffered = true;
          one.Buffer   = Asked(words);
        }
        else if (first == "offset")
        {
          std::string distance;
          std::string corner;
          words >> distance >> corner;
          one.Offset = {.Distance = std::strtod(distance.c_str(), nullptr), .Gaps = Gapped(corner)};
        }
        else if (first == "tolerance")
        {
          std::string length;
          words >> length;
          one.Within = {.Length = std::strtod(length.c_str(), nullptr), .Radians = 1e-10};
        }
        else if (first == "piece")
        {
          one.Answer.push_back(Drawn(words));
        }
        else
        {
          cases.push_back(one);
        }
      }
      return cases;
    }
  } // namespace

  TEST_CASE("reference: every offset is the reference program's, to the bit")
  {
    const std::vector<Reference> cases = References();
    REQUIRE(cases.size() == 200);
    for (const Reference &one : cases)
    {
      CAPTURE(one.Name);
      const auto got = one.Buffered ? Buffer(one.Shape, one.Buffer, one.Within)
                                    : Offset(one.Shape, one.Offset, one.Within);
      REQUIRE(got);
      CHECK(*got == one.Answer);
    }
  }

  TEST_CASE("reference: the same offset near the origin, at 5e6 and at 1e8")
  {
    for (const Reference &one : References())
    {
      CAPTURE(one.Name);
      const double distance = one.Buffered ? one.Buffer.Distance : one.Offset.Distance;
      const auto   found    = Offset(one.Shape, {.Distance = distance, .Gaps = Gap::Fillet},
                                     Tolerance::Of(Extents(one.Shape)));
      REQUIRE(found);
      for (const Vector by : {Vector{.X = 5e6, .Y = -5e6}, Vector{.X = -1e8, .Y = 1e8}})
      {
        Polyline there = one.Shape;
        for (Vertex &vertex : there.Vertices)
        {
          vertex.Position = vertex.Position + by;
        }
        const auto moved = Offset(there, {.Distance = distance, .Gaps = Gap::Fillet},
                                  Tolerance::Of(Extents(there)));
        REQUIRE(moved);
        REQUIRE(moved->size() == found->size());
        for (std::size_t index = 0; index < found->size(); ++index)
        {
          CHECK(std::abs(Length((*moved)[index]) - Length((*found)[index])) < 1e-6);
        }
      }
    }
  }
} // namespace Geoxy
