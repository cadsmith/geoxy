/***************************************************************************************************
 *  @file      source/met/solve.cpp
 *  @brief     Two parts met by the rule for their kinds, two runs, a run and a piece, two pieces,
 *             and any two parts by the rule for theirs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/solve.hpp"

#include "basis.hpp"
#include "kernel/piece.hpp"
#include "met/ellipse.hpp"
#include "met/parts.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "trig.hpp"

#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy::Met
{
  auto Evaluated(const Piece &piece, double along) -> Point
  {
    return Kernel::Sampled(piece, along).At;
  }

  auto Angled(const Matrix &unit, Point point) -> double
  {
    const Point seen = unit * point;
    return Trig::Atan2(seen.Y, seen.X);
  }

  namespace
  {
    /**
     *  How far a point is from a run's line, within about eight steps of a double at the size of
     *  the numbers: asked only whether it is within the tolerance, which is never finer than 32
     *  such steps, so an exact zero on the line would answer no differently.
     *
     *  @param  run   The run.
     *  @param  point The point.
     *  @return       The distance.
     */
    auto Astray(const Run &run, Point point) -> double
    {
      return std::abs(Cross(point - run.Start, run.Way)) / Length(run.Way);
    }

    /**
     *  Where two runs' lines cross, if it lies on both, or past either's end by the tolerance.
     *
     *  @param  first     One run, of some length.
     *  @param  second    The other.
     *  @param  tolerance How far past an end a crossing may lie.
     *  @return           The crossing; nothing where the lines are parallel or it lies off a run.
     */
    auto Athwart(const Run &first, const Run &second, const Tolerance &tolerance) -> Found
    {
      Found        found;
      const Vector way   = first.Way;
      const Vector other = second.Way;
      const double cross = Cross(way, other);
      if (Zero(cross))
      {
        return found;
      }
      const Vector gap    = second.Start - first.Start;
      const double along  = Cross(gap, other) / cross;
      const double beyond = Cross(gap, way) / cross;
      const double slack  = tolerance.Length / std::sqrt(Dot(way, way));
      const double spare  = tolerance.Length / Length(other);
      if (first.Along.Start - slack     <= along
          && along                      <= first.Along.End + slack
          && second.Along.Start - spare <= beyond
          && beyond                     <= second.Along.End + spare)
      {
        found.Locals.push_back({.First = along, .Second = beyond});
      }
      return found;
    }

  } // namespace

  auto Aligned(const Run &first, const Run &second, const Tolerance &tolerance)
      -> std::optional<Range>
  {
    const Vector way    = first.Way;
    const Vector other  = second.Way;
    const double square = Dot(way, way);
    /**
     *  A way shorter than the square root of the least double squares to nothing: the run is
     *  shorter than any tolerance the curves' numbers allow, and lies along nothing over a stretch
     *  longer than it.
     */
    if (!(square > 0.0))
    {
      return std::nullopt;
    }
    const double base  = Dot(second.Start - first.Start, way) / square;
    const double rate  = Dot(other, way) / square;
    const double from  = base + (rate * second.Along.Start);
    const double until = base + (rate * second.Along.End);
    const Range  span  = {.Start = std::max(first.Along.Start, std::min(from, until)),
                          .End   = std::min(first.Along.End, std::max(from, until))};
    if (!((span.End - span.Start) * Length(way) > tolerance.Length))
    {
      return std::nullopt;
    }
    const bool endless = !std::isfinite(span.Start) || !std::isfinite(span.End);
    if (endless)
    {
      const bool parallel =
          Orient(Point{}, Point{} + way, Point{} + other) == Orientation::Collinear;
      return parallel && Astray(first, second.Start) <= tolerance.Length
               ? std::optional<Range>{span}
               : std::nullopt;
    }
    const bool near = Astray(second, first.Start + (way * span.Start)) <= tolerance.Length
                   && Astray(second, first.Start + (way * span.End))   <= tolerance.Length;
    return near ? std::optional<Range>{span} : std::nullopt;
  }

  auto Against(const Run &first, const Run &second, const Tolerance &tolerance) -> Found
  {
    Found        found;
    const Vector way    = first.Way;
    const Vector other  = second.Way;
    const double square = Dot(way, way);
    if (Zero(square) || Zero(Dot(other, other)))
    {
      return found;
    }
    if (const std::optional<Range> span = Aligned(first, second, tolerance))
    {
      const double base = Dot(second.Start - first.Start, way) / square;
      const double rate = Dot(other, way) / square;
      const double low  = span->Start;
      const double high = span->End;
      const auto   back = [base, rate](double along) -> double
      {
        return (along - base) / rate;
      };
      const double start = back(low);
      const double end   = back(high);
      /**
       *  The second's points at the stretch's ends, whose feet on the first's line its ends are,
       *  lie farther from the first's line than the first's from the second's by the secant of
       *  their angle: past twice the tolerance, the first's within it, the runs part at more than
       *  sixty degrees, and the stretch is but the second's shadow, as where two runs cross nearly
       *  square and the first's short stretch about the crossing lies near the second's line, the
       *  second running far off. They cross there, if anywhere.
       */
      const bool shadow = Astray(first, second.Start + (other * start)) > 2.0 * tolerance.Length
                       || Astray(first, second.Start + (other * end))   > 2.0 * tolerance.Length;
      if (shadow)
      {
        return Athwart(first, second, tolerance);
      }
      found.Shared.push_back(
          {.First  = {.Start = low, .End = high},
           .Second = {.Start = std::min(start, end), .End = std::max(start, end)},
           .Along  = rate > 0.0});
      return found;
    }
    return Athwart(first, second, tolerance);
  }

  namespace
  {
    /**
     *  A stretch of a piece: the piece whole, and the part of its parameter meant.
     */
    struct Stretch
    {
      Piece Whole = {}; /**< The piece.                    */
      Range Part  = {}; /**< The stretch of its parameter. */
    };

    /**
     *  A point and how near to it.
     */
    struct Target
    {
      Point  At     = {};  /**< The point. */
      double Within = 0.0; /**< How near.  */
    };

    /**
     *  A place on a piece within a length of a point, sought by Newton's method on the squared
     *  distance from a guess: any such place shows the piece comes that near.
     *
     *  @param  piece  The piece.
     *  @param  target The point and the length.
     *  @param  guess  Where to start, from 0 to 1.
     *  @return        The place; nothing where a few steps do not find one.
     */
    auto Reached(const Piece &piece, Target target, double guess) -> std::optional<double>
    {
      double along = guess;
      for (int step = 0; step < 6; ++step)
      {
        const Kernel::Sample sample = Kernel::Sampled(piece, along);
        const Vector         apart  = sample.At - target.At;
        if (Length(apart) <= target.Within)
        {
          return along;
        }
        const double speed = Dot(sample.Way, sample.Way);
        if (!(speed > 0.0))
        {
          return std::nullopt;
        }
        along = std::clamp(along - (Dot(apart, sample.Way) / speed), 0.0, 1.0);
      }
      return std::nullopt;
    }

    /**
     *  Whether a stretch of a piece lies along another piece: at nine places of it, the other's
     *  nearest point within the tolerance, the middle tried first, then the ends, as a stretch that
     *  leaves the other mostly does so there. The places are taken on the piece whole, never on a
     *  part split from it, whose control points have lost bits far from the origin. Whether what
     *  lies along is a stretch shared or the band around a touch is told by where it ends. A place
     *  of the other within the tolerance, found by a few steps of Newton's method from the foot of
     *  the nearest place done, settles a place; the nearest point is sought only where they fail.
     *
     *  @param  stretch   The stretch.
     *  @param  other     The other piece.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Along(const Stretch &stretch, const Piece &other, const Tolerance &tolerance) -> bool
    {
      constexpr int                  Places = 8;
      constexpr std::array<int, 9>   Order  = {4, 0, 8, 2, 6, 1, 3, 5, 7};
      std::array<double, Places + 1> feet   = {};
      std::array<bool, Places + 1>   seen   = {};
      for (const int place : Order)
      {
        const auto   index = static_cast<std::size_t>(place);
        const double along = Kernel::Across(stretch.Part, static_cast<double>(place) / Places);
        const Point  point = Evaluated(stretch.Whole, along);
        std::optional<double> foot;
        for (std::size_t gap = 1; gap <= Places && !foot; ++gap)
        {
          for (const std::size_t near : {index - gap, index + gap})
          {
            if (!foot && near <= Places && seen.at(near))
            {
              foot = Reached(other, {.At = point, .Within = tolerance.Length}, feet.at(near));
            }
          }
        }
        if (!foot)
        {
          const double closest = Kernel::Closest(other, point);
          if (!(Distance(point, Evaluated(other, closest)) <= tolerance.Length))
          {
            return false;
          }
          foot = closest;
        }
        feet.at(index) = *foot;
        seen.at(index) = true;
      }
      return true;
    }

    /**
     *  Whether a piece is straight to within an angle: every step between two of its control points
     *  within the angle of its chord and running its way, so that the piece, in their hull, runs
     *  along it from end to end and never back; a piece that ends where it starts is not. Asked of
     *  whole pieces only, whose control points are the curve's own: a part split from one far from
     *  the origin has its control points on the grid of doubles there, and may look straight where
     *  the piece is not.
     *
     *  @param  piece   The piece.
     *  @param  radians The angle.
     *  @return         Whether it is.
     */
    auto Steady(const Piece &piece, double radians) -> bool
    {
      const auto   poles = Kernel::Poles(piece);
      const Vector way   = Projected(poles.back()) - Projected(poles.front());
      if (Zero(Dot(way, way)))
      {
        return false;
      }
      for (std::size_t index = 1; index < poles.size(); ++index)
      {
        const Vector step = Projected(poles[index]) - Projected(poles[index - 1]);
        if (std::abs(Cross(step, way)) > Trig::Sin(radians) * Length(step) * Length(way)
            || Dot(step, way)          < 0.0)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  The stretches a piece lying on a run's line shares with it, one for each pass of the piece
     *  along the line between the places it turns back, as a piece may run past its ends and back:
     *  each pass is where it runs one way only, so the place on it for a place on the run is found
     *  by halving within the pass, never on the piece whole, whose other passes cover it too.
     *
     *  @param run   The run, of some length.
     *  @param piece The piece, on the run's line.
     *  @param found What was found so far, the stretches added.
     */
    void Passes(const Run &run, const Piece &piece, Found &found)
    {
      const double square = Dot(run.Way, run.Way);
      const auto   level  = [&run, &piece, square](double along) -> double
      {
        return Dot(Evaluated(piece, along) - run.Start, run.Way) / square;
      };
      List<double> cuts =
          Kernel::Zeroed(piece, {.Through = run.Start, .Way = {.X = -run.Way.Y, .Y = run.Way.X}},
                         Kernel::Budget)
              .Turns;
      cuts.insert(cuts.begin(), 0.0);
      cuts.push_back(1.0);
      std::ranges::sort(cuts);
      for (std::size_t index = 1; index < cuts.size(); ++index)
      {
        const Range  pass  = {.Start = cuts[index - 1], .End = cuts[index]};
        const double from  = level(pass.Start);
        const double until = level(pass.End);
        const double low   = std::max(run.Along.Start, std::min(from, until));
        const double high  = std::min(run.Along.End, std::max(from, until));
        if (!(low < high))
        {
          continue;
        }
        const auto back = [&pass, &level, from, until](double target) -> double
        {
          Range span = pass;
          for (int step = 0; step < 64; ++step)
          {
            const double middle              = span.Start + ((span.End - span.Start) / 2);
            const bool   before              = (level(middle) < target) == (from < until);
            (before ? span.Start : span.End) = middle;
          }
          return span.Start + ((span.End - span.Start) / 2);
        };
        const bool   rising = from < until;
        const double head   = rising ? pass.Start : pass.End;
        const double tail   = rising ? pass.End : pass.Start;
        const double start  = Exactly(low, std::min(from, until)) ? head : back(low);
        const double end    = Exactly(high, std::max(from, until)) ? tail : back(high);
        found.Shared.push_back(
            {.First  = {.Start = low, .End = high},
             .Second = {.Start = std::min(start, end), .End = std::max(start, end)},
             .Along  = from < until});
      }
    }

    /**
     *  Whether two pieces part at more than an angle all along, one side of the other: every step
     *  between two control points of one turns from every step of the other by more than the angle,
     *  the same way round. A rational piece with positive weights runs at every parameter along a
     *  positive sum of its steps, and the cross product of two such sums outweighs the angle's
     *  share of their lengths where every pair of steps does; so the two cross at most once, a
     *  chord through two crossings running along both, and every crossing passes the same way.
     *
     *  @param  first   One piece.
     *  @param  second  The other.
     *  @param  radians The angle.
     *  @return         Whether they part.
     */
    auto Parting(const Piece &first, const Piece &second, double radians) -> bool
    {
      const double               sine   = Trig::Sin(radians);
      const auto                 poles  = Kernel::Poles(first);
      const auto                 others = Kernel::Poles(second);
      std::array<std::size_t, 2> sides  = {};
      for (std::size_t index = 1; index < poles.size(); ++index)
      {
        const Vector way = Projected(poles[index]) - Projected(poles[index - 1]);
        for (std::size_t other = 1; other < others.size(); ++other)
        {
          const Vector back   = Projected(others[other]) - Projected(others[other - 1]);
          const double turn   = Cross(back, way);
          const double apart  = sine * Length(way) * Length(back);
          sides[0]           += static_cast<std::size_t>(turn > apart);
          sides[1]           += static_cast<std::size_t>(turn < -apart);
        }
      }
      const std::size_t all = (poles.size() - 1) * (others.size() - 1);
      return sides[0] == all || sides[1] == all;
    }

    /**
     *  Where a piece crosses the lines a length either side of a run's, with its ends, in order and
     *  each once: between two of them it lies wholly within that length of the run's line or wholly
     *  outside.
     *
     *  @param  run    The run.
     *  @param  piece  The piece.
     *  @param  across How far either side, square to the run.
     *  @return        The places, from 0 to 1.
     */
    auto Marks(const Run &run, const Piece &piece, Vector across) -> List<double>
    {
      List<double> marks = {0.0, 1.0};
      for (const Vector shift : {across, -across})
      {
        const Kernel::Zeros zeros =
            Kernel::Zeroed(piece, {.Through = run.Start + shift, .Way = run.Way}, Kernel::Budget);
        marks.insert(marks.end(), zeros.Crossing.Roots.begin(), zeros.Crossing.Roots.end());
      }
      std::ranges::sort(marks);
      marks.erase(std::ranges::unique(marks).begin(), marks.end());
      return marks;
    }

    /**
     *  How far along a run a stretch of a piece reaches, each way: from its ends and the places it
     *  turns back along the run within it, the turns of where it crosses the line square to the
     *  run, so that the reach is exact, a cusp's tip and all.
     *
     *  @param  run     The run.
     *  @param  piece   The piece.
     *  @param  stretch The stretch.
     *  @return         The least place along the run and the most.
     */
    auto Reach(const Run &run, const Piece &piece, Range stretch) -> std::array<double, 2>
    {
      const double square = Dot(run.Way, run.Way);
      const auto   place  = [&run, square](Point point) -> double
      {
        return Dot(point - run.Start, run.Way) / square;
      };
      const double start  = place(Evaluated(piece, stretch.Start));
      const double end    = place(Evaluated(piece, stretch.End));
      double       least  = std::min(start, end);
      double       most   = std::max(start, end);
      const Vector across = {.X = -run.Way.Y, .Y = run.Way.X};
      for (const double turn :
           Kernel::Zeroed(piece, {.Through = run.Start, .Way = across}, Kernel::Budget).Turns)
      {
        if (stretch.Start < turn && turn < stretch.End)
        {
          least = std::min(least, place(Evaluated(piece, turn)));
          most  = std::max(most, place(Evaluated(piece, turn)));
        }
      }
      return {least, most};
    }

    /**
     *  The stretches of a piece that lie along a run's line, as two pieces' would be found: where
     *  the piece stays between the lines the tolerance either side, worked out from where it
     *  crosses them, more than four tolerances long, and not parting from the line at more than the
     *  angle, as a curve all but turned back at a cusp does not; each kept only where the run
     *  reaches all of it. A band among them, ending where the two part, is one hit, its kind told
     *  by whether they are tangent in it, whichever way either curve is spelled.
     *
     *  @param  run       The run.
     *  @param  piece     The piece.
     *  @param  tolerance The tolerance.
     *  @return           The stretches, the run first.
     */
    auto Lain(const Run &run, const Piece &piece, const Tolerance &tolerance) -> List<Overlap>
    {
      List<Overlap> lain;
      const double  square = Dot(run.Way, run.Way);
      const Vector  across =
          Vector{.X = -run.Way.Y, .Y = run.Way.X} * (tolerance.Length / std::sqrt(square));
      const List<double> marks = Marks(run, piece, across);
      const auto         off   = [&run, square](Point point) -> double
      {
        return std::abs(Cross(point - run.Start, run.Way)) / std::sqrt(square);
      };
      const auto place = [&run, square](Point point) -> double
      {
        return Dot(point - run.Start, run.Way) / square;
      };
      const Piece line = Kernel::Straight(run.Start, run.Start + run.Way);

      for (std::size_t index = 0; index + 1 < marks.size(); ++index)
      {
        const double low = marks[index];
        if (!(off(Evaluated(piece, low + ((marks[index + 1] - low) / 2))) <= tolerance.Length))
        {
          continue;
        }
        const double end         = marks[index + 1];
        const auto [least, most] = Reach(run, piece, {.Start = low, .End = end});
        const double length      = std::sqrt(square);
        const double slack       = tolerance.Length / length;
        const Piece  part        = Kernel::Split(Kernel::Split(piece, end)[0], low / end)[1];
        if (!((most - least) * length > 4.0 * tolerance.Length)
            || Parting(line, part, tolerance.Radians)
            || least < run.Along.Start - slack
            || most > run.Along.End + slack)
        {
          continue;
        }
        lain.push_back({.First  = {.Start = std::max(least, run.Along.Start),
                                   .End   = std::min(most, run.Along.End)},
                        .Second = {.Start = low, .End = end},
                        .Along  = place(Evaluated(piece, end)) > place(Evaluated(piece, low))});
      }
      return lain;
    }

  } // namespace

  auto Against(const Run &run, const Piece &piece, const Tolerance &tolerance) -> Found
  {
    Found        found;
    const double square = Dot(run.Way, run.Way);
    if (Zero(square))
    {
      return found;
    }
    const Run                  chord = {.Start = Kernel::Start(piece),
                                        .Way   = Kernel::End(piece) - Kernel::Start(piece),
                                        .Along = {.Start = 0, .End = 1}};
    const std::optional<Range> lying =
        Steady(piece, tolerance.Radians) ? Aligned(run, chord, tolerance) : std::nullopt;
    if (lying)
    {
      const double from  = Kernel::Closest(piece, run.Start + (run.Way * lying->Start));
      const double until = Kernel::Closest(piece, run.Start + (run.Way * lying->End));
      found.Shared.push_back(
          {.First  = *lying,
           .Second = {.Start = std::min(from, until), .End = std::max(from, until)},
           .Along  = from < until});
      return found;
    }
    const Kernel::Zeros zeros =
        Kernel::Zeroed(piece, {.Through = run.Start, .Way = run.Way}, Kernel::Budget);
    if (zeros.On)
    {
      Passes(run, piece, found);
      return found;
    }
    for (const Overlap &overlap : Lain(run, piece, tolerance))
    {
      found.Shared.push_back(overlap);
    }
    const auto add = [&found, &run, &piece, square](double along) -> void
    {
      const Point point = Evaluated(piece, along);
      found.Locals.push_back({.First = Dot(point - run.Start, run.Way) / square, .Second = along});
    };
    for (const double along : zeros.Crossing.Roots)
    {
      add(along);
    }
    for (const double along : zeros.Turns)
    {
      add(along);
    }
    return found;
  }

  namespace
  {
    /**
     *  How far a point is from a line: to the foot on it, or to the nearer end.
     *
     *  @param  point The point.
     *  @param  line  The line.
     *  @return       The distance to its nearest point; to its start where it has no length.
     */
    auto Clearance(Point point, const Line &line) -> double
    {
      const Vector way    = line.End - line.Start;
      const double square = Dot(way, way);
      const double along =
          square > 0.0 ? std::clamp(Dot(point - line.Start, way) / square, 0.0, 1.0) : 0.0;
      return Distance(point, Lerp(line.Start, line.End, along));
    }

    /**
     *  A piece's chord, and how far the piece strays from it.
     */
    struct Chord
    {
      Line   Ends  = {}; /**< The line from the piece's start to its end.                      */
      double Stray = 0;  /**< How far its farthest control point, and so it, is from the line. */
    };

    /**
     *  A piece's chord: the piece lies in its control points' hull, and every point within a
     *  distance of a line makes a hull that is, so it strays no farther than they do.
     *
     *  @param  piece The piece.
     *  @return       Its chord.
     */
    auto Chorded(const Piece &piece) -> Chord
    {
      Chord chord = {.Ends  = {.Start = Kernel::Start(piece), .End = Kernel::End(piece)},
                     .Stray = 0};
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        chord.Stray = std::max(chord.Stray, Clearance(Projected(pole), chord.Ends));
      }
      return chord;
    }

  } // namespace

  auto Apart(const Line &one, const Line &other) -> double
  {
    const auto side = [](const Line &line, Point point) -> int
    {
      return static_cast<int>(Orient(line.Start, line.End, point));
    };
    const std::array<int, 4> sides = {side(one, other.Start), side(one, other.End),
                                      side(other, one.Start), side(other, one.End)};
    const bool collinear = std::ranges::all_of(sides, [](int turn) -> bool { return turn == 0; });
    if (!collinear && sides[0] * sides[1] <= 0 && sides[2] * sides[3] <= 0)
    {
      return 0.0;
    }
    return std::min({Clearance(one.Start, other), Clearance(one.End, other),
                     Clearance(other.Start, one), Clearance(other.End, one)});
  }

  namespace
  {
    /**
     *  How far a line strays from another: the farther of its ends, as a point running along it is
     *  nowhere farther from the other than at one of them.
     *
     *  @param  one   The line.
     *  @param  other The other.
     *  @return       The distance.
     */
    auto Beyond(const Line &one, const Line &other) -> double
    {
      return std::max(Clearance(one.Start, other), Clearance(one.End, other));
    }

    /**
     *  What subdivision does next with two pieces.
     */
    enum class Choice : std::uint8_t
    {
      Drop,        /**< They cannot meet: their chords are too far apart.      */
      Intersect,   /**< Their chords' crossing is a candidate, polished later. */
      HalveFirst,  /**< The first is the larger: halve it.                     */
      HalveSecond, /**< The second is: halve it.                               */
    };

    /**
     *  How subdivision weighs two pieces.
     */
    struct Rule
    {
      double Length  = 0.0;   /**< The tolerance's length.                         */
      bool   Parting = false; /**< Whether the pieces part at more than its angle. */
    };

    /**
     *  What to do with two pieces that do not share a stretch: drop them where their chords are
     *  farther apart than the tolerance and how far the pieces stray from them; meet their chords
     *  where both are small and within the tolerance, where they part and so cross at most once, or
     *  where both are as straight as a double tells; else halve the larger.
     *
     *  @param  first  One piece.
     *  @param  second The other.
     *  @param  rule   The tolerance's length, and whether they part.
     *  @return        The step.
     */
    auto Stepped(const Piece &first, const Piece &second, Rule rule) -> Choice
    {
      const Chord  one   = Chorded(first);
      const Chord  two   = Chorded(second);
      const double stray = one.Stray + two.Stray;
      if (Apart(one.Ends, two.Ends) > rule.Length + stray)
      {
        return Choice::Drop;
      }
      const double reach     = std::max(Beyond(one.Ends, two.Ends), Beyond(two.Ends, one.Ends));
      const bool   together  = reach + stray <= rule.Length;
      const bool   single    = one.Stray <= rule.Length && two.Stray <= rule.Length && rule.Parting;
      const Vector firstSize = Size(Kernel::Hull(first));
      const Vector secondSize = Size(Kernel::Hull(second));
      const double firstSpan  = std::max(firstSize.X, firstSize.Y);
      const double secondSpan = std::max(secondSize.X, secondSize.Y);
      const double fine       = 1e-7 * std::max({firstSpan, secondSpan, 1e-300});
      if (together || single || (Kernel::Flat(first, fine) && Kernel::Flat(second, fine)))
      {
        return Choice::Intersect;
      }
      return firstSpan >= secondSpan ? Choice::HalveFirst : Choice::HalveSecond;
    }

    /**
     *  Where on a piece a point of its chord is: the parameter that part of the chord, where the
     *  piece's point there lies within a sixteenth of the chord of the point; else the place on the
     *  piece nearest the point. A piece whose weights differ much runs fast where they are light
     *  and slow where they are heavy, that part of it a long way from that part of its chord, too
     *  far for Newton's method to polish a candidate from.
     *
     *  @param  piece The piece.
     *  @param  along The part of its chord.
     *  @return       The parameter, from 0 to 1.
     */
    auto Reckoned(const Piece &piece, double along) -> double
    {
      const Point  start = Kernel::Start(piece);
      const Vector chord = Kernel::End(piece) - start;
      const Point  point = start + (chord * along);
      if (Distance(Kernel::Reached(piece, along), point) <= Length(chord) / 16.0)
      {
        return along;
      }
      return Kernel::Closest(piece, point);
    }

    /**
     *  Where two pieces' chords cross, as a parameter of each, clamped to the chord and taken to
     *  the piece; halfway along both where they are parallel.
     *
     *  @param  first  One piece.
     *  @param  second The other.
     *  @return        The two parameters.
     */
    auto Crossed(const Piece &first, const Piece &second) -> std::array<double, 2>
    {
      const Point  start    = Kernel::Start(first);
      const Point  other    = Kernel::Start(second);
      const Vector way      = Kernel::End(first) - start;
      const Vector beside   = Kernel::End(second) - other;
      const double crossing = Cross(way, beside);
      if (Zero(crossing))
      {
        return {0.5, 0.5};
      }
      return {Reckoned(first, std::clamp(Cross(other - start, beside) / crossing, 0.0, 1.0)),
              Reckoned(second, std::clamp(Cross(other - start, way) / crossing, 0.0, 1.0))};
    }

    /**
     *  A part of each of two pieces, as subdivision holds them.
     */
    struct Pair
    {
      Piece First;       /**< A part of the first piece.                      */
      Range FirstAlong;  /**< Its stretch of the first's.                     */
      Piece Second;      /**< A part of the second piece.                     */
      Range SecondAlong; /**< Its stretch of the second's.                    */
      bool  Off;         /**< Whether the part was found not to lie along it. */
    };

    /**
     *  The stretch two parts of straight pieces share, where their chords lie on one line within
     *  the tolerance over more than its length; its ends taken to the pieces whole.
     *
     *  @param  top       The parts.
     *  @param  wholes    The pieces whole, the first's and the second's.
     *  @param  tolerance The tolerance.
     *  @return           The stretch, or nothing.
     */
    auto Lined(const Pair &top, std::array<const Piece *, 2> wholes, const Tolerance &tolerance)
        -> std::optional<Overlap>
    {
      const Run chord                 = {.Start = Kernel::Start(top.First),
                                         .Way   = Kernel::End(top.First) - Kernel::Start(top.First),
                                         .Along = {.Start = 0, .End = 1}};
      const Run beside                = {.Start = Kernel::Start(top.Second),
                                         .Way   = Kernel::End(top.Second) - Kernel::Start(top.Second),
                                         .Along = {.Start = 0, .End = 1}};
      const std::optional<Range> span = Aligned(chord, beside, tolerance);
      if (!span)
      {
        return std::nullopt;
      }
      const Point  low   = chord.Start + (chord.Way * span->Start);
      const Point  high  = chord.Start + (chord.Way * span->End);
      const double from  = Kernel::Closest(*wholes[1], low);
      const double until = Kernel::Closest(*wholes[1], high);
      return Overlap{.First  = {.Start = Kernel::Closest(*wholes[0], low),
                                .End   = Kernel::Closest(*wholes[0], high)},
                     .Second = {.Start = std::min(from, until), .End = std::max(from, until)},
                     .Along  = from < until};
    }

    /**
     *  A pair of parts halved: the first part, when it is the larger, each half met with the second
     *  whole, tried afresh for lying along; else the second part, the first marked off where it was
     *  found not to lie along the other.
     *
     *  @param  top   The parts.
     *  @param  first Whether to halve the first.
     *  @param  tried Whether the first was tried for lying along the other and found not to.
     *  @return       The two pairs.
     */
    auto Halved(const Pair &top, bool first, bool tried) -> std::array<Pair, 2>
    {
      if (first)
      {
        const std::array<Piece, 2> halves = Kernel::Split(top.First, 0.5);
        const double               middle =
            top.FirstAlong.Start + ((top.FirstAlong.End - top.FirstAlong.Start) / 2);
        return {Pair{.First       = halves[0],
                     .FirstAlong  = {.Start = top.FirstAlong.Start, .End = middle},
                     .Second      = top.Second,
                     .SecondAlong = top.SecondAlong,
                     .Off         = false},
                Pair{.First       = halves[1],
                     .FirstAlong  = {.Start = middle, .End = top.FirstAlong.End},
                     .Second      = top.Second,
                     .SecondAlong = top.SecondAlong,
                     .Off         = false}};
      }
      const std::array<Piece, 2> halves = Kernel::Split(top.Second, 0.5);
      const double               middle =
          top.SecondAlong.Start + ((top.SecondAlong.End - top.SecondAlong.Start) / 2);
      return {Pair{.First       = top.First,
                   .FirstAlong  = top.FirstAlong,
                   .Second      = halves[0],
                   .SecondAlong = {.Start = top.SecondAlong.Start, .End = middle},
                   .Off         = top.Off || tried},
              Pair{.First       = top.First,
                   .FirstAlong  = top.FirstAlong,
                   .Second      = halves[1],
                   .SecondAlong = {.Start = middle, .End = top.SecondAlong.End},
                   .Off         = top.Off || tried}};
    }

  } // namespace

  auto Against(const Piece &first, const Piece &second, const Tolerance &tolerance) -> Found
  {
    Found       found;
    List<Pair>  pending  = {{.First       = first,
                             .FirstAlong  = {.Start = 0, .End = 1},
                             .Second      = second,
                             .SecondAlong = {.Start = 0, .End = 1},
                             .Off         = false}};
    std::size_t budget   = Kernel::Budget;
    const bool  straight = Steady(first, tolerance.Radians) && Steady(second, tolerance.Radians);
    while (!pending.empty())
    {
      const Pair top = pending.back();
      pending.pop_back();
      const Box  near   = Padded(Kernel::Hull(top.First), tolerance.Length);
      const Box  hull   = Kernel::Hull(top.Second);
      const bool shared = std::ranges::any_of(found.Shared,
                                              [&top](const Overlap &overlap) -> bool
                                              {
                                                return overlap.First.Start <= top.FirstAlong.Start
                                                    && top.FirstAlong.End  <= overlap.First.End;
                                              });
      if (!Overlaps(near, hull) || shared)
      {
        continue;
      }
      const double span    = Distance(Kernel::Start(top.First), Kernel::End(top.First));
      const bool   parting = Parting(top.First, top.Second, tolerance.Radians);
      const bool   tried   = span > 4.0 * tolerance.Length
                          && !top.Off
                          && !parting
                          && Contains(Padded(hull, tolerance.Length), Evaluated(top.First, 0.5));
      const bool lies = tried && Along({.Whole = first, .Part = top.FirstAlong}, second, tolerance);
      if (lies)
      {
        const double from  = Kernel::Closest(second, Kernel::Start(top.First));
        const double until = Kernel::Closest(second, Kernel::End(top.First));
        found.Shared.push_back(
            {.First  = top.FirstAlong,
             .Second = {.Start = std::min(from, until), .End = std::max(from, until)},
             .Along  = from < until});
        continue;
      }
      const Choice step =
          Stepped(top.First, top.Second, {.Length = tolerance.Length, .Parting = parting});
      if (step == Choice::Drop)
      {
        continue;
      }
      const std::optional<Overlap> lined =
          straight ? Lined(top, {&first, &second}, tolerance) : std::nullopt;
      if (lined)
      {
        found.Shared.push_back(*lined);
        continue;
      }
      if (step == Choice::Intersect)
      {
        const auto [along, beyond] = Crossed(top.First, top.Second);
        found.Locals.push_back(
            {.First = top.FirstAlong.Start + (along * (top.FirstAlong.End - top.FirstAlong.Start)),
             .Second =
                 top.SecondAlong.Start + (beyond * (top.SecondAlong.End - top.SecondAlong.Start))});
        continue;
      }
      if (budget == 0)
      {
        found.Spent = true;
        continue;
      }
      --budget;
      const std::array<Pair, 2> halves = Halved(top, step == Choice::HalveFirst, tried);
      pending.insert(pending.end(), halves.begin(), halves.end());
    }
    return found;
  }

  auto Matched(const Part &first, const Part &second, const Tolerance &tolerance) -> Found
  {
    return std::visit(
        [&tolerance](const auto &one, const auto &other) -> Found
        {
          if constexpr (requires { Against(one, other, tolerance); })
          {
            return Against(one, other, tolerance);
          }
          else
          {
            return Swapped(Against(other, one, tolerance));
          }
        },
        first.Shape, second.Shape);
  }
} // namespace Geoxy::Met
