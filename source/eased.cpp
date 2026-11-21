/***************************************************************************************************
 *  @file      source/eased.cpp
 *  @brief     A spline cut finer where its weights make it stiff, to be met in its place.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "eased.hpp"

#include "basis.hpp"
#include "kernel/piece.hpp"
#include "twofold.hpp"

#include "geoxy/spline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How many times a piece is halved at most: past every double's reach.
     */
    constexpr int DeepestHalving = 1070;

    /**
     *  A piece of a spline's span on its way to being cut: the piece, where it starts within the
     *  span from 0 to 1, how much of the span it is, and how often it was halved.
     */
    struct Cutting
    {
      Kernel::Piece     Shape = {};  /**< The piece.                               */
      Kernel::Hodograph Speed = {};  /**< Its derivative's numerator.              */
      Twofold::Double   Start = {};  /**< How far past the span's start it starts. */
      Twofold::Double   Rest  = {};  /**< How far short of the span's end it ends. */
      double            Width = 1.0; /**< How much of the span it is.              */
      int               Depth = 0;   /**< How many halvings made it.               */
    };

    /**
     *  Whether a piece's weights are all within so many times the lighter of its end weights.
     *
     *  @param  piece The piece.
     *  @param  mild  How many times.
     *  @return       Whether they are.
     */
    auto Balanced(const Kernel::Piece &piece, double mild) -> bool
    {
      const std::span<const Homogeneous> poles = std::span{piece.Poles}.first(piece.Degree + 1);
      const double                       light = std::min(poles.front().W, poles.back().W);
      return std::ranges::all_of(poles, [light, mild](const Homogeneous &pole) -> bool
                                 { return pole.W <= mild * light; });
    }

    /**
     *  A span's piece cut into even ones, in order: halved, the halves halved again, until each is
     *  even.
     *
     *  @param  whole The span's piece, and its derivative's numerator where the pieces are to carry
     *                theirs.
     *  @param  mild  Within how many times the lighter of its end weights an even piece's are.
     *  @return       The pieces, and where each lies within the span; none where halving as often
     *                as doubles allow leaves a piece uneven, its weights too far apart for doubles
     *                to hold beside each other.
     */
    auto Evened(const Kernel::Tread &whole, double mild) -> std::vector<Cutting>
    {
      std::vector<Cutting> even;
      std::vector<Cutting> pending = {Cutting{.Shape = whole.Shape,
                                              .Speed = whole.Speed,
                                              .Start = {},
                                              .Rest  = {},
                                              .Width = 1.0,
                                              .Depth = 0}};
      while (!pending.empty())
      {
        const Cutting top = pending.back();
        pending.pop_back();
        if (Balanced(top.Shape, mild))
        {
          even.push_back(top);
          continue;
        }
        if (top.Depth == DeepestHalving)
        {
          return {};
        }
        const std::array<Kernel::Piece, 2>     halves = Kernel::Split(top.Shape, 0.5);
        const std::array<Kernel::Hodograph, 2> speeds = Kernel::Halved(top.Speed);
        const double                           half   = top.Width / 2.0;
        pending.push_back({.Shape = halves[1],
                           .Speed = speeds[1],
                           .Start = Twofold::Add(top.Start, Twofold::Lifted(half)),
                           .Rest  = top.Rest,
                           .Width = half,
                           .Depth = top.Depth + 1});
        pending.push_back({.Shape = halves[0],
                           .Speed = speeds[0],
                           .Start = top.Start,
                           .Rest  = Twofold::Add(top.Rest, Twofold::Lifted(half)),
                           .Width = half,
                           .Depth = top.Depth + 1});
      }
      return even;
    }
  } // namespace

  auto Eased(const Spline &spline, double mild) -> std::optional<Easing>
  {
    if (spline.Weights.empty() || !Valid(spline))
    {
      return std::nullopt;
    }
    const auto [lightest, heaviest] = std::ranges::minmax(spline.Weights);
    if (heaviest <= mild * lightest)
    {
      return std::nullopt;
    }
    const auto degree = static_cast<std::size_t>(spline.Degree);
    Easing easing = {.Shape  = {.Degree = spline.Degree, .Knots = {}, .Poles = {}, .Weights = {}},
                     .Pieces = {}};
    bool   cut    = false;
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (!(spline.Knots[span] < spline.Knots[span + 1]))
      {
        continue;
      }
      const int                  power = Leveling(spline, span);
      const std::vector<Cutting> pieces =
          Evened({.Shape = Kernel::Section(spline, span), .Speed = {}}, mild);
      if (pieces.empty())
      {
        return Easing{};
      }
      cut = cut || pieces.size() > 1;
      for (const Cutting &piece : pieces)
      {
        easing.Pieces.push_back({.Span  = span,
                                 .From  = spline.Knots[span],
                                 .To    = spline.Knots[span + 1],
                                 .Start = piece.Start,
                                 .Rest  = piece.Rest,
                                 .Width = piece.Width});
        for (std::size_t index = easing.Shape.Poles.empty() ? 0 : 1; index <= degree; ++index)
        {
          const Homogeneous &pole = piece.Shape.Poles.at(index);
          easing.Shape.Poles.push_back({.X = pole.X / pole.W, .Y = pole.Y / pole.W});
          easing.Shape.Weights.push_back(std::ldexp(pole.W, -power));
        }
      }
    }
    if (!cut)
    {
      return std::nullopt;
    }
    const std::size_t count = easing.Pieces.size();
    for (std::size_t knot = 0; knot <= count; ++knot)
    {
      const std::size_t times = knot == 0 || knot == count ? degree + 1 : degree;
      easing.Shape.Knots.insert(easing.Shape.Knots.end(), times, static_cast<double>(knot));
    }
    return easing;
  }

  auto Back(const Easing &easing, Twofold::Double along) -> Sited
  {
    const auto            last  = static_cast<double>(easing.Pieces.size() - 1);
    const double          floor = std::clamp(std::floor(along.High), 0.0, last);
    const Slice          &piece = easing.Pieces[static_cast<std::size_t>(floor)];
    const Twofold::Double part  = Twofold::Difference(along, Twofold::Lifted(floor));
    const Twofold::Double left  = Twofold::Difference(Twofold::Lifted(1.0), part);
    const Twofold::Double width = Twofold::Lifted(piece.Width);
    const Twofold::Double span =
        Twofold::Difference(Twofold::Lifted(piece.To), Twofold::Lifted(piece.From));
    return {
        .Span  = piece.Span,
        .Ahead = Twofold::Multiply(Twofold::Add(piece.Start, Twofold::Multiply(part, width)), span),
        .Short = Twofold::Multiply(Twofold::Add(piece.Rest, Twofold::Multiply(left, width)), span)};
  }

  auto Paved(const Spline &spline) -> std::optional<Paving>
  {
    if (spline.Weights.empty() || !Valid(spline))
    {
      return std::nullopt;
    }
    const auto [lightest, heaviest] = std::ranges::minmax(spline.Weights);
    if (heaviest <= Gentle * lightest)
    {
      return std::nullopt;
    }
    const auto degree = static_cast<std::size_t>(spline.Degree);
    Paving     paving = {};
    bool       cut    = false;
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (!(spline.Knots[span] < spline.Knots[span + 1]))
      {
        continue;
      }
      const std::vector<Cutting> pieces = Evened(
          {.Shape = Kernel::Section(spline, span), .Speed = Kernel::Derived(spline, span)}, Gentle);
      if (pieces.empty())
      {
        return std::nullopt;
      }
      cut = cut || pieces.size() > 1;
      for (const Cutting &piece : pieces)
      {
        paving.Treads.push_back({.Shape = piece.Shape, .Speed = piece.Speed});
        paving.Pieces.push_back({.Span  = span,
                                 .From  = spline.Knots[span],
                                 .To    = spline.Knots[span + 1],
                                 .Start = piece.Start,
                                 .Rest  = piece.Rest,
                                 .Width = piece.Width});
      }
    }
    return cut ? std::optional<Paving>{std::move(paving)} : std::nullopt;
  }

  auto Homed(const Slice &piece, double along) -> double
  {
    const Twofold::Double within = Twofold::Add(
        piece.Start, Twofold::Multiply(Twofold::Lifted(along), Twofold::Lifted(piece.Width)));
    const Twofold::Double span =
        Twofold::Difference(Twofold::Lifted(piece.To), Twofold::Lifted(piece.From));
    const Twofold::Double place =
        Twofold::Add(Twofold::Lifted(piece.From), Twofold::Multiply(within, span));
    return place.High + place.Low;
  }

  auto Homed(const Easing &easing, double along) -> double
  {
    const Sited           site  = Back(easing, Twofold::Lifted(along));
    const auto            last  = static_cast<double>(easing.Pieces.size() - 1);
    const double          floor = std::clamp(std::floor(along), 0.0, last);
    const Slice          &piece = easing.Pieces[static_cast<std::size_t>(floor)];
    const Twofold::Double place = Twofold::Add(Twofold::Lifted(piece.From), site.Ahead);
    return place.High + place.Low;
  }
} // namespace Geoxy
