/***************************************************************************************************
 *  @file      source/graze.cpp
 *  @brief     Two ready curves that cross shallowly, nearly touch or nearly miss, a line or an arc
 *             against an arc, told from what they keep as Intersect tells them, to the bit, or left
 *             to Intersect.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "graze.hpp"

#include "abut.hpp"
#include "axis.hpp"
#include "chord.hpp"
#include "crossing.hpp"
#include "elliptic.hpp"
#include "lean.hpp"
#include "met/linear.hpp"
#include "real.hpp"
#include "tangency.hpp"
#include "thorough.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  The gap, in tolerances, below which Intersect's band takes two curves for touching, and
     *  twice the depth below which it takes two crossing for one touch: an eighth short of each
     *  bound.
     */
    constexpr double Below = 7.0 / 8.0;

    /**
     *  The gap, in tolerances, past which Intersect's band takes two curves for apart, and twice
     *  the depth past which it takes two crossing for two crossings: an eighth past each bound.
     */
    constexpr double Atop = 9.0 / 8.0;

    /**
     *  How far past a line's ends a point on its carrier may lie and be its, as Intersect places a
     *  crossing: a few steps of a double.
     */
    constexpr double Overrun = 16.0 * std::numeric_limits<double>::epsilon();

    /**
     *  How far apart two crossings' places along the first curve must lie, as fractions of a line
     *  or quarter turns round an arc, to be put in Intersect's order by what is worked out here:
     *  far past what Intersect's roundings of either place can take.
     */
    constexpr double Spread = 0x1p-36;

    /**
     *  The angle, as a sine, past which a direction is clear of an arc's ends for Intersect's angle
     *  tests: per radian of the arc's start and sweep, and per radius of its center from the
     *  origin, far past what the library's sine, cosine and arc tangent and the sums of angles
     *  leave.
     */
    constexpr double Angular = 0x1p-40;

    /**
     *  The square of the sine of Steep's angle's complement, a little short: a carrier whose
     *  center's foot lies within it times the radius cuts the circle at a sine past Steep by more
     *  than Roundly's roundings take, with no root taken.
     */
    constexpr double Squat = (1.0 - (Met::Steep * Met::Steep)) * (1.0 - 0x1p-20);

    /**
     *  What Intersect's band tells of two curves near each other.
     */
    enum class Tidings : std::uint8_t
    {
      Unsure,    /**< Nothing told here: Intersect decides. */
      Nothing,   /**< They meet in nothing.                 */
      Touch,     /**< They touch once, at the point given.  */
      Crossings, /**< They cross at the rough points given. */
      Steep,     /**< They cross at a steep angle, an end near the other curve, as Roundly keeps
                      them: crossings at the rough points given, at an end where an end kept lies
                      within the tolerance of one. */
      Shy        /**< Roundly gives them up, an end near the other curve, for Shy to tell. */
    };

    /**
     *  What Intersect's band tells of two curves, and what it is told from.
     */
    struct Band
    {
      Tidings              Kind    = Tidings::Unsure; /**< What it tells.                */
      Point                At      = {};              /**< The touch, where it is one.   */
      std::array<Point, 2> Rough   = {};              /**< The crossings, roughly.       */
      std::size_t          Count   = 0;               /**< How many crossings there are. */
      bool                 Parting = false;           /**< Whether Intersect first asks for the
                                                           carriers' touch, and goes on only
                                                           where they certainly do not.    */
      int Lost = 0;                                   /**< How many of Intersect's tests on the way
                                                           are not told here: where any is
                                                           not, Intersect decides.              */
      double               Length  = 0.0;             /**< The tolerance. */
      std::array<Point, 2> Carried = {};              /**< A line's carrier's two crossings with the
                                                           circle, roughly, the nearer its start first. */
      std::array<double, 2> Places = {};              /**< Their places along the line. */
      std::array<Point, 4>  Ends   = {};              /**< The ends kept at a crossing. */
      std::size_t           Ended  = 0;               /**< How many there are.          */
      bool                  Near   = false; /**< Whether an end lies within four tolerances
                                                 of the other curve, as Bend asks.          */
      std::array<double, 2> Asked = {};     /**< Where Intersect asks the curves to touch: their
                                                 places, roughly, as it gives them.         */

      /**
       *  One more crossing, roughly.
       *
       *  @param point Its point.
       */
      void Add(Point point)
      {
        Rough.at(Count) = point;
        ++Count;
      }
    };

    /**
     *  A ready line and a ready arc, as Intersect's direct way reads them: the line's way from its
     *  start to its end and its square, each rounded as plain arithmetic gives it.
     */
    struct Bow
    {
      const Glancing::Held *Line   = nullptr; /**< The line.                                  */
      const Glancing::Held *Round  = nullptr; /**< The arc.                                   */
      Vector                Way    = {};    /**< The line's way.                                  */
      double                Square = 0.0;   /**< Its square.                                      */
      double                Slack  = 0.0;   /**< The slack at their size, as Leeway gives it.     */
      bool                  Steep  = false; /**< Whether a steep crossing is told, or only bands. */
      const Glancing::Held *First =
          nullptr; /**< Which of the two is first, as Intersect takes them. */
    };

    /**
     *  Where a ready arc's center's foot lies on a ready line's carrier, as Roundly works it out:
     *  its place along the line, and its distance from the center.
     */
    struct Plumb
    {
      double Middle = 0.0; /**< The place.    */
      double Off    = 0.0; /**< The distance. */
    };

    /**
     *  A length Intersect compares distances with, and a slack past what its roundings of a
     *  distance can take at the curves' size: a distance past the length and the slack is past it
     *  as Intersect works it out, one short of the length less the slack short of it.
     */
    struct Ambit
    {
      double Length = 0.0; /**< The length. */
      double Slack  = 0.0; /**< The slack.  */
    };

    /**
     *  The slack of a Reach at a size: Intersect works a distance between points of the curves out,
     *  their difference, its length, a foot on a line and a radius taken from it, within 16 u S of
     *  the true one, S the largest magnitude of their numbers; the slack is eight times that.
     *
     *  @param  size S.
     *  @return      The slack.
     */
    auto Leeway(double size) -> double
    {
      return 0x1p-46 * size;
    }

    /**
     *  The arc a ready arc is.
     *
     *  @param  round The ready arc.
     *  @return       Its arc.
     */
    auto Arced(const Glancing::Held &round) -> const Arc &
    {
      return std::get<Arc>(*round.Shape);
    }

    /**
     *  How Intersect's band tells two curves by the gap where they come nearest, in tolerances,
     *  less than zero where they cross, as Gapped tells them.
     *
     *  @param  apart The gap.
     *  @return       Nothing past Above; a touch short of Below, or crossing no deeper than twice
     *                Below; crossings deeper than twice Above; unsure within an eighth of a bound.
     */
    auto Gapped(double apart) -> Tidings
    {
      const double deep = -apart / 2.0;
      if (apart > Atop)
      {
        return Tidings::Nothing;
      }
      if (apart >= Below || (deep >= Below && deep <= Atop))
      {
        return Tidings::Unsure;
      }
      return deep < Below ? Tidings::Touch : Tidings::Crossings;
    }
  } // namespace

  struct Glancing::Bands
  {
    /**
     *  Whether a ready arc is a whole circle.
     *
     *  @param  round The ready arc.
     *  @return       Whether it is.
     */
    static auto Whole(const Held &round) -> bool { return std::abs(round.Sweep) >= 2.0 * Pi; }

    /**
     *  How far a point is from a ready arc's circle, as Intersect's Off works it out.
     *
     *  @param  point The point.
     *  @param  round The ready arc.
     *  @return       The distance.
     */
    static auto Off(Point point, const Held &round) -> double
    {
      const Arc &arc = Arced(round);
      return std::abs(Distance(point, arc.Center) - arc.Radius);
    }

    /**
     *  How far a point is from a ready line, as Intersect's Away works it out for a line of some
     *  length: to its foot, or to its nearer end.
     *
     *  @pre          The line's square is above zero.
     *  @param  point The point.
     *  @param  bow   The line, and the arc.
     *  @return       The distance.
     */
    static auto Away(Point point, const Bow &bow) -> double
    {
      const Point  start = bow.Line->Ends[0];
      const double along = std::clamp(Dot(point - start, bow.Way) / bow.Square, 0.0, 1.0);
      return Distance(point, start + (bow.Way * along));
    }

    /**
     *  Whether a point of a ready arc's circle lies within its sweep, as Intersect's Swept tells
     *  it, on the same ends by the same steps: by the turns from its ends seen from the center.
     *
     *  @param  round The ready arc.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    static auto Swept(const Held &round, Point point) -> bool
    {
      const double size = std::abs(round.Sweep);
      if (size >= 2.0 * Pi)
      {
        return true;
      }
      const Arc   &arc     = Arced(round);
      const bool   forward = round.Sweep > 0.0;
      const Vector start   = (forward ? round.Ends[0] : round.Ends[1]) - arc.Center;
      const Vector end     = (forward ? round.Ends[1] : round.Ends[0]) - arc.Center;
      const Vector toward  = point - arc.Center;
      if (size <= Pi)
      {
        return Cross(start, toward) >= 0.0 && Cross(toward, end) >= 0.0;
      }
      const bool rest = Cross(end, toward) > 0.0 && Cross(toward, start) > 0.0;
      return !rest;
    }

    /**
     *  Whether a way from a ready arc's center points within its sweep, as Intersect's angle tests
     *  tell it, its heading's turns from the start counted: 1 within, -1 off, 0 where they may
     *  answer either way.
     *
     *  The ends kept lie within a few steps of a double of the arc's circle at the start's and the
     *  end's angles, a few steps at the size of the center's numbers over the radius as an angle;
     *  the heading of the way, its turn from the start and the count of turns are each within a few
     *  steps at the size of the angles. Each turn from an end is told by its sign where its sine
     *  lies past Angular times all of those sizes, far past what they leave: both are needed to
     *  tell a way within a sweep up to a half turn, one to tell it off, and the other way round for
     *  a longer sweep. A whole circle's way is within it unless it points near its seam.
     *
     *  @param  round  The ready arc.
     *  @param  toward The way.
     *  @return        Which.
     */
    static auto Within(const Held &round, Vector toward) -> int
    {
      const Arc   &arc     = Arced(round);
      const double opening = arc.Start.HeldInDegrees()
                               ? std::abs(arc.Start.InDegrees()) * (Pi / 90.0)
                               : std::abs(arc.Start.InRadians());
      const double turns   = ((4.0 + opening + std::abs(round.Sweep)) * arc.Radius)
                           + std::abs(arc.Center.X)
                           + std::abs(arc.Center.Y);
      const double edge    = Angular * turns * (std::abs(toward.X) + std::abs(toward.Y));
      if (Whole(round))
      {
        const Vector seam = round.Ends[0] - arc.Center;
        return std::max(static_cast<int>(std::abs(Cross(seam, toward)) > edge),
                        static_cast<int>(Dot(seam, toward) < 0.0));
      }
      const auto   back  = static_cast<std::size_t>(!(round.Sweep > 0.0));
      const Vector start = round.Ends.at(back) - arc.Center;
      const Vector end   = round.Ends.at(1 - back) - arc.Center;
      const double left  = Cross(start, toward);
      const double right = Cross(toward, end);
      const int    clear =
          static_cast<int>(std::abs(left) > edge) * static_cast<int>(std::abs(right) > edge);
      const int narrow  = static_cast<int>(std::abs(round.Sweep) <= Pi);
      const int both    = static_cast<int>(left > 0.0) * static_cast<int>(right > 0.0);
      const int neither = static_cast<int>(left < 0.0) * static_cast<int>(right < 0.0);
      const int inside  = (narrow * both) + ((1 - narrow) * (1 - neither));
      return clear * ((2 * inside) - 1);
    }

    /**
     *  How far round a ready arc a point of it lies, from its start the way it sweeps: a number
     *  that grows with the place along it, no angle taken.
     *
     *  @param  round The ready arc.
     *  @param  point The point.
     *  @return       From 0 at the start to 4 a whole turn on.
     */
    static auto Turned(const Held &round, Point point) -> double
    {
      const Arc   &arc    = Arced(round);
      const Vector start  = round.Ends[0] - arc.Center;
      const Vector toward = point - arc.Center;
      const double side   = round.Sweep > 0.0 ? Cross(start, toward) : Cross(toward, start);
      const double level  = Dot(start, toward);
      const double share  = level / (std::abs(level) + std::abs(side));
      return side >= 0.0 ? 1.0 - share : 3.0 + share;
    }

    /**
     *  Where Intersect places a crossing's point on a ready line or arc, told: 1 where it places
     *  it, a line's by the fraction of its way it works out, here by the same steps, within a few
     *  steps of a double of its ends, an arc's where Within tells its angle within the sweep; -1
     *  where it does not; 0 where Intersect's own place is to be asked. With it a number in the
     *  order of the place Intersect gives: a line's fraction, an arc's turn from its start.
     *
     *  @param  curve The ready curve.
     *  @param  point The point.
     *  @return       Whether, and its place's order.
     */
    static auto Placed(const Held &curve, Point point) -> std::pair<int, double>
    {
      if (const auto *line = std::get_if<Line>(&*curve.Shape))
      {
        const Vector way   = line->End - line->Start;
        const double along = Dot(point - line->Start, way) / Dot(way, way);
        const int    placed =
            static_cast<int>(-Overrun <= along) * static_cast<int>(along <= 1.0 + Overrun);
        return {(2 * placed) - 1, along};
      }
      return {Within(curve, point - Arced(curve).Center), Turned(curve, point)};
    }

    /**
     *  Whether a ready arc points a way from its center, as Intersect's Facing tells it: as Within
     *  tells it, or where Within does not tell it, by Facing's own steps.
     *
     *  @param  round  The ready arc.
     *  @param  toward The way.
     *  @return        Whether it does.
     */
    static auto Facing(const Held &round, Vector toward) -> bool
    {
      const int within = Within(round, toward);
      if (within != 0)
      {
        return within > 0;
      }
      return Earliest(Elliptic(Arced(round)), Heading(toward).InRadians()).has_value();
    }

    /**
     *  Of the carriers' two crossings, the one nearer a rough crossing, as Intersect picks it: by
     *  the two distances from the doubles nearest the crossings, rounded, M + H where they are
     *  equal, as its quick roots, its short sums and its exact way each pick it.
     *
     *  @param  ends The carriers' two crossings, M + H and M - H.
     *  @param  near The rough crossing.
     *  @return      The nearer.
     */
    static auto Nearer(const std::array<Chord::Settle, 2> &ends, Point near) -> Point
    {
      return Distance(ends[0].At, near) <= Distance(ends[1].At, near) ? ends[0].At : ends[1].At;
    }

    /**
     *  Crossings found as Intersect's Exactly finds them: each rough crossing's exact one, how many
     *  steps are not told, how many points Intersect places on neither curve or one, and their
     *  places' order along the first curve, Intersect's own where asked.
     */
    struct Found
    {
      Abutment              Hits   = {}; /**< The crossings, in the order found.   */
      std::array<double, 2> Places = {}; /**< Their places' order along the first. */
      int                   Lost   = 0;  /**< How many steps are not told: a point not
                                              certain.                                  */
      int Off = 0;                       /**< How many Intersect places on neither curve or
                                              one, so that its Exactly finds nothing.   */
    };

    /**
     *  Whether Intersect places a point on an arc whose ends At gives are held, told without
     *  placing it: where the point lies within eight steps of a double at its size, no less than 1
     *  nor the radius, of either end. Crossing's Located places a point on an arc by its angle,
     *  else by its nearest place, kept where that lies inside the arc or within that length of the
     *  point; and the nearest place is weighed against both ends, as At gives them, so lies no
     *  farther from the point than either, the distances worked out alike either way round.
     *
     *  @param  curve The curve, held.
     *  @param  point The point.
     *  @return       Whether it is told so; where not, Intersect's own place is asked.
     */
    static auto Tipped(const Held &curve, Point point) -> bool
    {
      if (!curve.Tipped)
      {
        return false;
      }
      const double size =
          std::max({1.0, std::abs(point.X), std::abs(point.Y), std::abs(Arced(curve).Radius)});
      const double reach = 8.0 * size * std::numeric_limits<double>::epsilon();
      return Distance(curve.Tips[0], point) <= reach || Distance(curve.Tips[1], point) <= reach;
    }

    /**
     *  The crossings of a band made exact as Intersect makes them: each rough crossing the nearer
     *  of the carriers' two exact ones, as Nearer picks it, each certain; placed on both curves, as
     *  Placed tells it or, where it does not, as Intersect places them; their places in order along
     *  the first, Intersect's own where two lie within Spread of each other.
     *
     *  @param  curves The two ready curves.
     *  @param  band   The rough crossings.
     *  @param  ends   The carriers' two crossings, M + H and M - H.
     *  @return        The crossings found.
     */
    static auto Finding(const Holds &curves, const Band &band,
                        const std::array<Chord::Settle, 2> &ends) -> Found
    {
      Found found = {.Lost = static_cast<int>(!ends[0].Sure) + static_cast<int>(!ends[1].Sure)};
      int   asked = 0;
      for (std::size_t index = 0; index < band.Count; ++index)
      {
        const Point point          = Nearer(ends, band.Rough.at(index));
        const auto [first, place]  = Placed(curves[0], point);
        const int second           = Placed(curves[1], point).first;
        found.Off                 += static_cast<int>(first < 0) + static_cast<int>(second < 0);
        asked                     += static_cast<int>(first == 0) + static_cast<int>(second == 0);
        found.Places.at(found.Hits.Count)     = place;
        found.Hits.Marks.at(found.Hits.Count) = {.At = point, .Kind = Contact::Cross};
        ++found.Hits.Count;
      }
      /**
       *  Where an angle is not told, or two places are near each other, Intersect's own places are
       *  asked, and the hits put in its order by them: its sort keeps two at one place as found.
       */
      const int close  = static_cast<int>(found.Hits.Count == 2)
                       * static_cast<int>(!(std::abs(found.Places[1] - found.Places[0]) > Spread));
      const int tipped = static_cast<int>(found.Hits.Count == 1)
                       * static_cast<int>(Tipped(curves[0], found.Hits.Marks[0].At))
                       * static_cast<int>(Tipped(curves[1], found.Hits.Marks[0].At));
      if (asked + close > 0 && tipped == 0)
      {
        const std::array<const Curve *, 2> shapes = {curves[0].get().Shape, curves[1].get().Shape};
        for (std::size_t index = 0; index < found.Hits.Count; ++index)
        {
          const std::optional<Crossings::Crossed> placed =
              Crossings::Placed(shapes, {}, found.Hits.Marks.at(index).At);
          found.Off              += static_cast<int>(!placed.has_value());
          found.Places.at(index)  = placed.value_or(Crossings::Crossed{}).Along[0];
        }
      }
      if (found.Hits.Count == 2 && found.Places[1] < found.Places[0])
      {
        std::swap(found.Hits.Marks[0], found.Hits.Marks[1]);
        std::swap(found.Places[0], found.Places[1]);
      }
      return found;
    }

    /**
     *  The crossings of a band made exact as Intersect makes them, as Finding finds them, in order
     *  along the first curve.
     *
     *  @param  curves The two ready curves.
     *  @param  band   The rough crossings.
     *  @param  ends   The carriers' two crossings, M + H and M - H.
     *  @return        The hits; nothing where any step is not told, or Intersect places a point off
     *                 either curve.
     */
    static auto Settled(const Holds &curves, const Band &band,
                        const std::array<Chord::Settle, 2> &ends) -> std::optional<Abutment>
    {
      const Found found = Finding(curves, band, ends);
      if (found.Lost + found.Off > 0)
      {
        return std::nullopt;
      }
      return found.Hits;
    }

    /**
     *  Whether a point lies farther than a length from a ready arc's circle, as Intersect's Off
     *  tells it: by the square of its distance from the center, rounded, where that lies past the
     *  radius and the length and the slack, or short of the radius less them, or within the length
     *  less the slack of the radius, the slack far past what the square's roundings take; else by
     *  Off.
     *
     *  @param  point The point.
     *  @param  round The ready arc.
     *  @param  clear The length, and the slack.
     *  @return       Whether it does.
     */
    static auto Beyond(Point point, const Held &round, Ambit clear) -> bool
    {
      const Arc   &arc    = Arced(round);
      const Vector apart  = point - arc.Center;
      const double square = Dot(apart, apart);
      const double inner  = arc.Radius - (clear.Length + clear.Slack);
      const double outer  = arc.Radius + (clear.Length + clear.Slack);
      const double near   = arc.Radius - (clear.Length - clear.Slack);
      const double far    = arc.Radius + (clear.Length - clear.Slack);
      if (square < inner * inner * Lean::Tight || square > outer * outer * Lean::Loose)
      {
        return true;
      }
      if (square > near * near * Lean::Loose && square < far * far * Lean::Tight)
      {
        return false;
      }
      return Off(point, round) > clear.Length;
    }

    /**
     *  Whether a point lies farther than a length from a ready line, as Intersect's Away tells it:
     *  where its distance from the line's carrier, as a cross product with the way squared,
     *  rounded, lies past the length and the slack, the slack far past what that product's
     *  roundings take, as the distance from the line is no less; else by Away.
     *
     *  @param  point The point.
     *  @param  bow   The line, and the arc.
     *  @param  clear The length, and the slack.
     *  @return       Whether it does.
     */
    static auto Distant(Point point, const Bow &bow, Ambit clear) -> bool
    {
      const double across = Cross(point - bow.Line->Ends[0], bow.Way);
      const double far    = clear.Length + clear.Slack;
      if (across * across > far * far * bow.Square * Lean::Loose)
      {
        return true;
      }
      return Away(point, bow) > clear.Length;
    }

    /**
     *  The line's ends each farther than a length from the arc's circle, as Intersect's Off tells
     *  it.
     *
     *  @param  bow   The line and the arc.
     *  @param  clear The length, and the slack.
     *  @return       Whether they are.
     */
    static auto Ended(const Bow &bow, Ambit clear) -> bool
    {
      return Beyond(bow.Line->Ends[0], *bow.Round, clear)
          && Beyond(bow.Line->Ends[1], *bow.Round, clear);
    }

    /**
     *  The arc's ends each farther than a length from the line, as Intersect's Away tells it: both
     *  kept ends where all are asked, else only an arc's, a whole circle's seam being no end.
     *
     *  @param  bow   The line and the arc.
     *  @param  clear The length, and the slack.
     *  @param  all   Whether a whole circle's seam is asked too.
     *  @return       Whether they are.
     */
    static auto Cleared(const Bow &bow, Ambit clear, bool all) -> bool
    {
      const Held &round = *bow.Round;
      if (!all && Whole(round))
      {
        return true;
      }
      return Distant(round.Ends[0], bow, clear) && Distant(round.Ends[1], bow, clear);
    }

    /**
     *  Each arc's ends farther than a length from the other's circle, as Intersect's Off tells it:
     *  both kept ends where all are asked, else only an arc's, a whole circle's seam being no end.
     *
     *  @param  rounds The two ready arcs.
     *  @param  clear  The length, and the slack.
     *  @param  all    Whether a whole circle's seam is asked too.
     *  @return        Whether they are.
     */
    static auto Cleared(const Holds &rounds, Ambit clear, bool all) -> bool
    {
      for (std::size_t side = 0; side < 2; ++side)
      {
        const Held &round   = rounds.at(side);
        const Held &against = rounds.at(1 - side);
        if (!all && Whole(round))
        {
          continue;
        }
        if (!Beyond(round.Ends[0], against, clear) || !Beyond(round.Ends[1], against, clear))
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  The slack two ready arcs are told within, at the largest magnitude of their numbers.
     *
     *  @param  rounds The two ready arcs.
     *  @return        The slack.
     */
    static auto Leeway(const Holds &rounds) -> double
    {
      double size = 0.0;
      for (const Held &round : rounds)
      {
        const Arc &arc = Arced(round);
        size           = std::max(
            {size, std::abs(arc.Center.X) + arc.Radius, std::abs(arc.Center.Y) + arc.Radius});
      }
      return Geoxy::Leeway(size);
    }

    /**
     *  A ready line's carrier and a ready arc's circle as a quadratic in the place t along the
     *  line, P + t d against |X - O| = r, as Preparing's Chorded works it out: a = d.d, b = d.(P -
     *  O), and the discriminant b^2 - a c for c = |P - O|^2 - r^2, each in lean double-double
     *  within its bound; and how far past zero the discriminant must lie for Intersect's quick
     *  roots to find the carriers certainly apart.
     */
    struct Quadric
    {
      Lean::Double A     = {};  /**< a.                                               */
      Lean::Double B     = {};  /**< b.                                               */
      Lean::Double Disc  = {};  /**< b^2 - a c.                                       */
      double       OffA  = 0.0; /**< How far a may be off.                            */
      double       OffB  = 0.0; /**< How far b may be off.                            */
      double       OffD  = 0.0; /**< How far the discriminant may be off.             */
      double       Quick = 0.0; /**< How far past zero Intersect's must certainly be. */
    };

    /**
     *  A ready line's carrier and a ready arc's circle as a quadratic, worked out as Preparing's
     *  Chorded works it out, each bound its.
     *
     *  Intersect's quick roots work the discriminant out in double-double from the same exact
     *  differences, and take the carriers for certainly apart where it lies past four times its
     *  bound, twice the sum of 2^-98 (b^2 + |a c|), 2^-98 times 2 |b| B, a (2 G + r^2) and |c| a, B
     *  and G the size of b's and of |P - O|^2's terms, at most 2^-97 (3 B^2 + a (4 G + 3 r^2)) with
     *  |b| at most B and |c| at most G + r^2. Its value is within that bound of the true one, so a
     *  true discriminant past five times it, here past 2^-94 times that, leaves Intersect's past
     *  four times its own.
     *
     *  @param  bow The line and the arc.
     *  @return     The quadratic.
     */
    static auto Quadratic(const Bow &bow) -> Quadric
    {
      const Held        &line    = *bow.Line;
      const Held        &round   = *bow.Round;
      const Lean::Double a       = {.High = line.Square, .Low = line.Under};
      const Lean::Way    way     = {Lean::Double{.High = bow.Way.X, .Low = line.Slip.X},
                                    Lean::Double{.High = bow.Way.Y, .Low = line.Slip.Y}};
      const Lean::Way    gap     = Lean::Apart(line.Ends[0], Arced(round).Center);
      const Lean::Sized  b       = Lean::Dot(way, gap);
      const Lean::Sized  away    = Lean::Dot(gap, gap);
      const Lean::Double squared = {.High = round.Square, .Low = round.Under};
      const Lean::Double c       = Lean::Plus(away.Value, Lean::Minus(squared));
      const Lean::Double lead    = Lean::Times(b.Value, b.Value);
      const Lean::Double rest    = Lean::Times(a, c);
      const double       offA    = Lean::Grain * a.High;
      const double       offB    = Lean::Grain * b.Size;
      const double       offC    = Lean::Grain * (away.Size + squared.High);
      return {.A     = a,
              .B     = b.Value,
              .Disc  = Lean::Plus(lead, Lean::Minus(rest)),
              .OffA  = offA,
              .OffB  = offB,
              .OffD  = (((2.0 * std::abs(b.Value.High) * offB)
                         + (offB * offB)
                         + (a.High * offC)
                         + (std::abs(c.High) * offA)
                         + (offA * offC))
                        * Lean::Loose)
                     + (Lean::Grain * (lead.High + std::abs(rest.High))),
              .Quick = 0x1p-94
                     * ((3.0 * b.Size * b.Size)
                        + (a.High * ((4.0 * away.Size) + (3.0 * squared.High))))};
    }

    /**
     *  The two crossings of a ready line's carrier and a ready arc's circle, M + H then M - H, M =
     *  P + d (-b / a) and H = d (sqrt(b^2 - a c) / a), each certified as Preparing's Cut certifies
     *  it, by the same bounds but the root's: the discriminant's bound over the root of the
     *  discriminant less its bound, which holds wherever the true one is not below zero, and is no
     *  bound at all, so certifies nothing, where the discriminant is not clear above zero; a
     *  coordinate on an axis made certain as Chord's Zeroed makes it.
     *
     *  @param  bow   The line and the arc.
     *  @param  shape Their quadratic.
     *  @return       The crossings, each certain or not.
     */
    static auto Halved(const Bow &bow, const Quadric &shape) -> std::array<Chord::Settle, 2>
    {
      const double       room    = std::max((shape.Disc.High * Lean::Tight) - shape.OffD, 0.0);
      const Held        &line    = *bow.Line;
      const Lean::Way    way     = {Lean::Double{.High = bow.Way.X, .Low = line.Slip.X},
                                    Lean::Double{.High = bow.Way.Y, .Low = line.Slip.Y}};
      const Lean::Rooted root    = Lean::Root(shape.Disc);
      const double       inverse = line.Inverse;
      const Lean::Double middle  = Lean::Over(Lean::Minus(shape.B), shape.A, inverse);
      const Lean::Double half    = Lean::Over(root.Value, shape.A, inverse);
      const double       offRoot =
          (shape.OffD / std::sqrt(room) * Lean::Loose) + (Lean::Grain * root.Value.High);
      const double offMiddle =
          ((shape.OffB + (std::abs(middle.High) * shape.OffA)) * inverse * Lean::Loose)
          + (Lean::Grain * std::abs(middle.High));
      const double  offHalf = ((offRoot + (half.High * shape.OffA)) * inverse * Lean::Loose)
                            + (Lean::Grain * half.High);
      const Point   start   = line.Ends[0];
      Chord::Halves halves;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Lean::Double step = way.at(axis);
        halves.Middle.at(axis)  = Lean::Moved(axis == 0 ? start.X : start.Y, step, middle);
        halves.Half.at(axis)    = Lean::Times(step, half);
        halves.Off.at(axis)     = 2.0 * std::abs(step.High) * (offMiddle + offHalf);
      }
      const Line &straight = std::get<Line>(*line.Shape);
      const Arc  &arc      = Arced(*bow.Round);
      const auto  zero     = [&straight, &arc, &shape](Point spot, double halfway) -> bool
      {
        return static_cast<int>(Axis::Fits(shape.A.High * halfway))
                 * static_cast<int>(Axis::On(straight, spot))
                 * static_cast<int>(Axis::On(arc, spot))
             > 0;
      };
      return Chord::Zeroed(Chord::Ends(halves), halves, zero);
    }

    /**
     *  Two ready arcs' circles, as Preparing's Circled works them out: their centers' gap g exact,
     *  |g|^2, the chord's k = |g|^2 + r1^2 - r2^2 and big = 4 |g|^2 r1^2 - k^2, each in lean
     *  double-double within its bound; and how far past zero big must lie for Intersect's quick
     *  roots to find the circles certainly apart.
     */
    struct Circles
    {
      Lean::Way    Gap   = {};  /**< g.                                               */
      Lean::Double Apart = {};  /**< |g|^2.                                           */
      Lean::Double Chord = {};  /**< k.                                               */
      Lean::Double Big   = {};  /**< 4 |g|^2 r1^2 - k^2.                              */
      double       OffG  = 0.0; /**< How far |g|^2 may be off.                        */
      double       OffK  = 0.0; /**< How far k may be off.                            */
      double       OffB  = 0.0; /**< How far big may be off.                          */
      double       Quick = 0.0; /**< How far past zero Intersect's must certainly be. */
    };

    /**
     *  Two ready arcs' circles, worked out as Preparing's Circled works them out, each bound its.
     *
     *  Intersect's quick roots work big out in double-double from the same exact differences, and
     *  take the circles for certainly apart where it lies past four times its bound, twice the sum
     *  of 2^-98 (4 G r1^2 + k^2), 2^-98 times 4 r1^2 G and 2 |k| (2 G + r1^2 + r2^2), G = |g|^2:
     *  with S = G + r1^2 + r2^2, which |k| is at most and 2 |g| r1 at most too, at most 7 2^-97
     *  S^2. Its value is within that bound of the true one, so a true big past five times it, here
     *  past 2^-91 S^2, leaves Intersect's past four times its own.
     *
     *  @param  rounds The two ready arcs.
     *  @return        Their circles.
     */
    static auto Circled(const Holds &rounds) -> Circles
    {
      const Held        &one   = rounds[0];
      const Held        &other = rounds[1];
      const Lean::Way    gap   = Lean::Apart(Arced(other).Center, Arced(one).Center);
      const Lean::Sized  apart = Lean::Dot(gap, gap);
      const Lean::Double inner = {.High = one.Square, .Low = one.Under};
      const Lean::Double outer = {.High = other.Square, .Low = other.Under};
      const Lean::Double chord = Lean::Plus(Lean::Plus(apart.Value, inner), Lean::Minus(outer));
      const Lean::Double four  = Lean::Times(apart.Value, inner);
      const Lean::Double lead  = {.High = 4.0 * four.High, .Low = 4.0 * four.Low};
      const Lean::Double tail  = Lean::Times(chord, chord);
      const double       offG  = Lean::Grain * apart.Size;
      const double       offK  = offG + (Lean::Grain * (apart.Size + inner.High + outer.High));
      const double       size  = apart.Size + inner.High + outer.High;
      return {.Gap   = gap,
              .Apart = apart.Value,
              .Chord = chord,
              .Big   = Lean::Plus(lead, Lean::Minus(tail)),
              .OffG  = offG,
              .OffK  = offK,
              .OffB =
                  (((4.0 * inner.High * offG) + (2.0 * std::abs(chord.High) * offK) + (offK * offK))
                   * Lean::Loose)
                  + (Lean::Grain * (lead.High + tail.High)),
              .Quick = 0x1p-91 * size * size};
    }

    /**
     *  The two crossings of two ready arcs' circles, M + H then M - H, M = O1 + g k / (2 |g|^2) and
     *  H = g' sqrt(big) / (2 |g|^2), g' a quarter turn of g, each certified as Preparing's Circled
     *  certifies it, by the same bounds but the root's, which is taken as Halved takes it for a
     *  line; a coordinate on an axis made certain as Chord's Zeroed makes it.
     *
     *  @param  rounds The two ready arcs.
     *  @param  circle Their circles.
     *  @return        The crossings, each certain or not.
     */
    static auto Halved(const Holds &rounds, const Circles &circle) -> std::array<Chord::Settle, 2>
    {
      const double       room    = std::max((circle.Big.High * Lean::Tight) - circle.OffB, 0.0);
      const Arc         &first   = Arced(rounds[0]);
      const Lean::Way   &gap     = circle.Gap;
      const Lean::Double twice   = {.High = 2.0 * circle.Apart.High, .Low = 2.0 * circle.Apart.Low};
      const double       inverse = 1.0 / twice.High;
      const Lean::Rooted root    = Lean::Root(circle.Big);
      const Lean::Double along   = Lean::Over(circle.Chord, twice, inverse);
      const Lean::Double across  = Lean::Over(root.Value, twice, inverse);
      const double       offRoot =
          (circle.OffB / std::sqrt(room) * Lean::Loose) + (Lean::Grain * root.Value.High);
      const double offAlong =
          ((circle.OffK + (std::abs(along.High) * 2.0 * circle.OffG)) * inverse * Lean::Loose)
          + (Lean::Grain * std::abs(along.High));
      const double offAcross =
          ((offRoot + (across.High * 2.0 * circle.OffG)) * inverse * Lean::Loose)
          + (Lean::Grain * across.High);
      const Chord::Halves halves = {
          .Middle = {Lean::Moved(first.Center.X, gap[0], along),
                     Lean::Moved(first.Center.Y, gap[1], along)},
          .Half   = {Lean::Minus(Lean::Times(gap[1], across)), Lean::Times(gap[0], across)},
          .Off    = {
              2.0 * ((std::abs(gap[0].High) * offAlong) + (std::abs(gap[1].High) * offAcross)),
              2.0 * ((std::abs(gap[1].High) * offAlong) + (std::abs(gap[0].High) * offAcross))}};
      const Arc &second = Arced(rounds[1]);
      const auto zero   = [&first, &second, &twice](Point spot, double middle) -> bool
      {
        return static_cast<int>(Axis::Fits(twice.High * middle))
                 * static_cast<int>(Axis::On(first, spot))
                 * static_cast<int>(Axis::On(second, spot))
             > 0;
      };
      return Chord::Zeroed(Chord::Ends(halves), halves, zero);
    }

    /**
     *  A ready line and a ready arc whose carriers come within four tolerances of touching, told as
     *  Roundly's Touching, then Bend, tell them: nothing where the foot of the center lies off
     *  either curve, every end far; else Bend's band about the foot, the carriers certainly apart,
     *  its touch where the arc points to the foot, as Intersect's angle test tells it, halfway
     *  between the foot and the circle's point toward it.
     *
     *  @param  bow    The line and the arc.
     *  @param  foot   Where the center's foot lies on the line's carrier.
     *  @param  length The tolerance.
     *  @return        What Intersect tells of them.
     */
    static auto Touching(const Bow &bow, const Plumb &foot, double length) -> Band
    {
      const Held  &round  = *bow.Round;
      const Arc   &arc    = Arced(round);
      const Ambit  span   = {.Length = 4.0 * length, .Slack = bow.Slack};
      const double radius = arc.Radius;
      const Point  spot   = bow.Line->Ends[0] + (bow.Way * foot.Middle);
      if (!Ended(bow, span) || !Cleared(bow, span, false))
      {
        return {};
      }
      const bool onBoth = 0.0 < foot.Middle && foot.Middle < 1.0 && Swept(round, spot);
      if (!onBoth)
      {
        return {.Kind = Tidings::Nothing};
      }
      const bool first = bow.Line == bow.First;
      Band       band  = {.Kind = Gapped((foot.Off - radius) / length), .Parting = true};
      band.Asked       = first ? std::array{foot.Middle, 0.0} : std::array{0.0, foot.Middle};
      if (band.Kind == Tidings::Touch)
      {
        band.Lost = static_cast<int>(!Facing(round, spot - arc.Center));
        band.At   = Crossings::Nearing(*bow.Line->Shape, arc);
      }
      if (band.Kind == Tidings::Crossings)
      {
        const double half =
            std::sqrt((radius * radius) - (foot.Off * foot.Off)) / std::sqrt(bow.Square);
        band.Add(bow.Line->Ends[0] + (bow.Way * (foot.Middle - half)));
        band.Add(bow.Line->Ends[0] + (bow.Way * (foot.Middle + half)));
      }
      return band;
    }

    /**
     *  A ready line and a ready arc whose carriers cross at a sine under Steep, deeper than four
     *  tolerances, told as Bend tells them, every end far: its crossings, each kept where it lies
     *  within the line and the sweep.
     *
     *  @param  bow    The line and the arc.
     *  @param  foot   Where the center's foot lies on the line's carrier.
     *  @param  length The tolerance.
     *  @return        What Intersect tells of them.
     */
    static auto Shallow(const Bow &bow, const Plumb &foot, double length) -> Band
    {
      const Held  &round  = *bow.Round;
      const double radius = Arced(round).Radius;
      const Ambit  span   = {.Length = 4.0 * length, .Slack = bow.Slack};
      if (!Ended(bow, span) || !Cleared(bow, span, true))
      {
        return {};
      }
      const double half =
          std::sqrt((radius * radius) - (foot.Off * foot.Off)) / std::sqrt(bow.Square);
      Band band = {.Kind = Gapped((foot.Off - radius) / length)};
      for (const double along : {foot.Middle - half, foot.Middle + half})
      {
        const Point point = bow.Line->Ends[0] + (bow.Way * along);
        if (0.0 < along && along < 1.0 && Swept(round, point))
        {
          band.Add(point);
        }
      }
      return band;
    }

    /**
     *  Whether a point lies within half the tolerance of either of a line's carrier's crossings
     *  with a circle, as Settle and Shy ask it.
     *
     *  @param  point The point.
     *  @param  band  The crossings, and the tolerance.
     *  @return       Whether it does.
     */
    static auto Kept(Point point, const Band &band) -> bool
    {
      const double half = band.Length / 2;
      return std::max(static_cast<int>(Distance(point, band.Carried[0]) <= half),
                      static_cast<int>(Distance(point, band.Carried[1]) <= half))
           > 0;
    }

    /**
     *  A ready line and a ready arc whose carriers cross at a sine past Steep, told as Roundly
     *  tells them where its exact passes leave them, an end near the other curve or a whole
     *  circle's seam near the line: nothing where both crossings lie off the line and its ends past
     *  eight tolerances; else the line's and the arc's ends, each farther than four tolerances from
     *  the other curve or within half a tolerance of a crossing, kept; and each crossing on the
     *  line, or at an end of it kept, and on the arc, or at an end of it kept. Where an end is
     *  neither, or a whole circle's seam near the line is at no crossing, Roundly gives them up,
     *  for Shy.
     *
     *  @param  bow    The line and the arc.
     *  @param  foot   Where the center's foot lies on the line's carrier.
     *  @param  length The tolerance.
     *  @return        What Roundly tells of them.
     */
    static auto Pierced(const Bow &bow, const Plumb &foot, double length) -> Band
    {
      const Held  &round  = *bow.Round;
      const double radius = Arced(round).Radius;
      const double clear  = 4.0 * length;
      const Ambit  span   = {.Length = clear, .Slack = bow.Slack};
      const double half =
          std::sqrt((radius * radius) - (foot.Off * foot.Off)) / std::sqrt(bow.Square);
      Band band = {.Kind = Tidings::Steep, .Length = length};
      for (std::size_t side = 0; side < 2; ++side)
      {
        const double sign     = (2.0 * static_cast<double>(side)) - 1.0;
        band.Places.at(side)  = foot.Middle + (sign * half);
        band.Carried.at(side) = bow.Line->Ends[0] + (bow.Way * band.Places.at(side));
      }
      const auto inside = [&band](std::size_t side) -> int
      {
        return static_cast<int>(0.0 < band.Places.at(side))
             * static_cast<int>(band.Places.at(side) < 1.0);
      };
      const std::array<Point, 2> &ends  = bow.Line->Ends;
      const std::array<Point, 2> &cuts  = round.Ends;
      const Ambit                 wide  = {.Length = 2.0 * clear, .Slack = bow.Slack};
      const int                   aloof = (1 - inside(0))
                                        * (1 - inside(1))
                                        * static_cast<int>(Beyond(ends[0], round, wide))
                                        * static_cast<int>(Beyond(ends[1], round, wide));
      const std::array<bool, 4>   far = {Beyond(ends[0], round, span), Beyond(ends[1], round, span),
                                         Distant(cuts[0], bow, span), Distant(cuts[1], bow, span)};
      band.Near                       = static_cast<int>(far[0])
                                          * static_cast<int>(far[1])
                                          * static_cast<int>(far[2])
                                          * static_cast<int>(far[3])
                                     == 0;
      if (aloof > 0)
      {
        return {.Kind = Tidings::Nothing};
      }
      /**
       *  Each end near the other curve at a crossing is kept; a whole circle's seam at one is no
       *  end, and not kept.
       */
      const bool                 whole  = Whole(round);
      const std::array<Point, 4> points = {ends[0], ends[1], cuts[0], cuts[1]};
      int                        given  = 0;
      for (std::size_t index = 0; index < 4; ++index)
      {
        const Point point         = points.at(index);
        const int   near          = static_cast<int>(!far.at(index));
        const int   kept          = near * static_cast<int>(Kept(point, band));
        const int   seam          = static_cast<int>(index >= 2) * static_cast<int>(whole);
        given                    += near - kept;
        band.Ends.at(band.Ended)  = point;
        band.Ended               += static_cast<std::size_t>(kept * (1 - seam));
      }
      if (given > 0)
      {
        band.Kind = Tidings::Shy;
        return band;
      }
      /**
       *  Each crossing on the line, or within half a tolerance of either of its ends, and on the
       *  arc, or within half a tolerance of an end of it not a whole circle's seam.
       */
      for (std::size_t side = 0; side < 2; ++side)
      {
        const Point  point = band.Carried.at(side);
        const double reach = length / 2;
        const auto   met   = [point, reach](Point end) -> int
        {
          return static_cast<int>(Distance(end, point) <= reach);
        };
        const int onLine = std::max({inside(side), met(ends[0]), met(ends[1])});
        const int onArc =
            std::max({static_cast<int>(Swept(round, point)),
                      (1 - static_cast<int>(whole)) * std::max(met(cuts[0]), met(cuts[1]))});
        if (onLine * onArc > 0)
        {
          band.Add(point);
        }
      }
      return band;
    }

    /**
     *  A ready line and a ready arc near each other, told as Intersect's direct way tells them,
     *  every test on its own rough numbers made on the same numbers by the same steps: Roundly's
     *  foot of the center on the line's carrier and its distance from the center; then
     *  - where the carrier comes within four tolerances of the circle, as Touching tells it;
     *  - where it passes farther outside, nothing, every end lying farther than four tolerances
     *    from the other curve, as Roundly asks;
     *  - where it cuts the circle at a sine under Steep, as Shallow tells it.
     *
     *  @param  bow    The line and the arc.
     *  @param  length The tolerance.
     *  @return        What Intersect tells of them.
     */
    static auto Bent(const Bow &bow, double length) -> Band
    {
      const Arc   &arc    = Arced(*bow.Round);
      const double clear  = 4.0 * length;
      const double radius = arc.Radius;
      if (!(bow.Square > 0.0) || !(radius > 16.0 * clear))
      {
        return {};
      }
      const Point  start  = bow.Line->Ends[0];
      const double middle = Dot(arc.Center - start, bow.Way) / bow.Square;
      const Point  spot   = start + (bow.Way * middle);
      const Vector toward = spot - arc.Center;
      /**
       *  Where only bands are asked, a crossing plainly past Steep is left at once, no root taken.
       */
      if (!bow.Steep && Dot(toward, toward) < radius * radius * Squat)
      {
        return {};
      }
      const Plumb  foot  = {.Middle = middle, .Off = Distance(spot, arc.Center)};
      const double ratio = foot.Off / radius;
      if (!(std::abs(foot.Off - radius) > clear))
      {
        return Touching(bow, foot, length);
      }
      if (foot.Off > radius)
      {
        const Ambit span = {.Length = clear, .Slack = bow.Slack};
        const int   near = static_cast<int>(!(foot.Off - radius > 2.0 * clear));
        return {.Kind = Tidings::Nothing,
                .Lost = near
                      * (static_cast<int>(!Ended(bow, span))
                         + static_cast<int>(!Cleared(bow, span, true)))};
      }
      if ((1.0 - (ratio * ratio)) > Met::Steep * Met::Steep)
      {
        return bow.Steep ? Pierced(bow, foot, length) : Band{};
      }
      return Shallow(bow, foot, length);
    }

    /**
     *  Two ready arcs whose circles come near each other, told as Intersect's Rim tells them: their
     *  centers' gap and the gap between the circles where they come nearest, outside or within, the
     *  circles' points on the line of centers there, and where they cross, roughly, from the first
     *  circle; a touch halfway between those points where each arc points to its own, as
     *  Intersect's angle test tells it.
     *
     *  @param  rounds The two ready arcs.
     *  @param  length The tolerance.
     *  @return        What Intersect tells of them.
     */
    static auto Rim(const Holds &rounds, double length) -> Band
    {
      const Held  &one     = rounds[0];
      const Held  &other   = rounds[1];
      const Arc   &first   = Arced(one);
      const Arc   &second  = Arced(other);
      const double clear   = 4.0 * length;
      const double apart   = Distance(first.Center, second.Center);
      const double outside = apart - (first.Radius + second.Radius);
      const double inside  = std::abs(first.Radius - second.Radius) - apart;
      const bool   kissing = std::abs(outside) < std::abs(inside);
      const double gap     = kissing ? outside : inside;
      const bool   shallow = gap >= -clear;
      const Ambit  span    = {.Length = clear, .Slack = Leeway(rounds)};
      const int    lost    = static_cast<int>(!(apart > clear))
                           + static_cast<int>(!(gap <= clear))
                           + static_cast<int>(!Cleared(rounds, span, !shallow));
      if (lost > 0)
      {
        return {};
      }
      const Vector toward = (second.Center - first.Center) / apart;
      const bool   larger = first.Radius >= second.Radius;
      const double inward = kissing || larger ? 1.0 : -1.0;
      const double onward = kissing ? -1.0 : inward;
      const Vector near   = toward * inward;
      const Vector far    = toward * onward;
      const Point  here   = first.Center + (near * first.Radius);
      const Point  there  = second.Center + (far * second.Radius);
      /**
       *  Near tangent, where they come nearest lies on both or Intersect finds nothing; that it
       *  lies off either, which Kissing has asked of one point between those two, is left to
       *  Intersect.
       */
      const int swept = static_cast<int>(Swept(one, here)) * static_cast<int>(Swept(other, there));
      Band band = {.Kind = Gapped(gap / length), .Lost = static_cast<int>(shallow) * (1 - swept)};
      if (band.Kind == Tidings::Touch)
      {
        band.Lost += static_cast<int>(!Facing(one, near)) + static_cast<int>(!Facing(other, far));
        band.At    = Crossings::Nearing({&first, &second}, {inward, onward});
      }
      if (band.Kind == Tidings::Crossings)
      {
        const double inner  = first.Radius;
        const double outer  = second.Radius;
        const double along  = ((apart * apart) + (inner * inner) - (outer * outer)) / (2.0 * apart);
        const double height = std::sqrt(std::max(0.0, (inner * inner) - (along * along)));
        const Point  base   = first.Center + (toward * along);
        const Vector rise   = Perpendicular(toward) * height;
        for (const Point point : {base + rise, base - rise})
        {
          if (shallow || (Swept(one, point) && Swept(other, point)))
          {
            band.Add(point);
          }
        }
      }
      return band;
    }

    /**
     *  Two ready arcs whose circles cross at a sine past Steep, told as Circularly tells them where
     *  Preparing's exact pass leaves them, an end near the other circle: the crossings, roughly,
     *  from the smaller circle, as Circularly works them out; each arc's ends, each farther than
     *  four tolerances from the other's circle, or within half a tolerance of a crossing and kept;
     *  and each crossing on both arcs, within each sweep or within half a tolerance of an end.
     *  Where an end is neither, or a whole circle's seam lies within four tolerances of the other
     *  circle, Circularly gives them up, and so is Intersect left to decide.
     *
     *  @param  rounds The two ready arcs, of circles crossing at a sine past Steep, more than four
     *                 tolerances from touching.
     *  @param  length The tolerance.
     *  @return        What Circularly tells of them.
     */
    static auto Riven(const Holds &rounds, double length) -> Band
    {
      const Held  &one    = rounds[0];
      const Held  &other  = rounds[1];
      const Arc   &first  = Arced(one);
      const Arc   &second = Arced(other);
      const double clear  = 4.0 * length;
      const double apart  = Distance(first.Center, second.Center);
      const bool   ahead  = first.Radius <= second.Radius;
      const Arc   &small  = ahead ? first : second;
      const Arc   &large  = ahead ? second : first;
      const double less   = small.Radius;
      const double more   = large.Radius;
      const Vector toward = (large.Center - small.Center) / apart;
      const double along  = ((apart * apart) + (less * less) - (more * more)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (less * less) - (along * along)));
      const Point  base   = small.Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      Band         band   = {.Kind = Tidings::Steep, .Length = length};
      band.Carried        = {base + rise, base - rise};
      /**
       *  Each end near the other circle at a crossing is kept; one near it elsewhere, or a whole
       *  circle's seam near it anywhere, leaves the pair to Intersect.
       */
      for (std::size_t side = 0; side < 2; ++side)
      {
        const Held &round   = rounds.at(side);
        const Held &against = rounds.at(1 - side);
        const int   whole   = static_cast<int>(Whole(round));
        for (const Point end : round.Ends)
        {
          const int near            = static_cast<int>(!(Off(end, against) > clear));
          const int kept            = near * (1 - whole) * static_cast<int>(Kept(end, band));
          band.Lost                += near - kept;
          band.Ends.at(band.Ended)  = end;
          band.Ended               += static_cast<std::size_t>(kept);
        }
      }
      for (const Point point : band.Carried)
      {
        const auto met = [point, length](const Held &round) -> bool
        {
          return Distance(round.Ends[0], point) <= length / 2
              || Distance(round.Ends[1], point) <= length / 2;
        };
        if ((Swept(one, point) || met(one)) && (Swept(other, point) || met(other)))
        {
          band.Add(point);
        }
      }
      return band;
    }

    /**
     *  Two ready arcs near each other, told as Intersect's direct way tells them, every test on its
     *  own rough numbers made on the same numbers by the same steps:
     *  - Circularly, where the circles keep more than eight tolerances clear of each other, apart
     *    or one within the other, finds nothing;
     *  - where they come within four of touching, asks Kissing for a touch on the line of centers:
     *    nothing where it lies off either arc, every end far; else, the circles certainly apart,
     *    Rim's band;
     *  - where they cross at a sine under Steep, Rim's band; past it, where steep crossings are
     *    asked, as Riven tells them;
     *  - where they keep four to eight clear, nothing where every end lies farther than four
     *    tolerances from the other circle.
     *
     *  @param  rounds The two ready arcs.
     *  @param  length The tolerance.
     *  @param  steep  Whether a steep crossing is told, or only bands.
     *  @return        What Intersect tells of them.
     */
    static auto Ringed(const Holds &rounds, double length, bool steep) -> Band
    {
      const Held  &one    = rounds[0];
      const Held  &other  = rounds[1];
      const Arc   &first  = Arced(one);
      const Arc   &second = Arced(other);
      const double clear  = 4.0 * length;
      const double inner  = first.Radius;
      const double outer  = second.Radius;
      if (!(inner > 16.0 * clear)
          || !(outer > 16.0 * clear)
          || (first.Center == second.Center && Exactly(inner, outer)))
      {
        return {};
      }
      const double apart   = Distance(first.Center, second.Center);
      const double sum     = inner + outer;
      const double gap     = std::abs(inner - outer);
      const bool   crossed = apart < sum && apart > gap;
      if (!crossed && std::min(std::abs(apart - sum), std::abs(apart - gap)) > 2.0 * clear)
      {
        return {.Kind = Tidings::Nothing};
      }
      if (!(std::abs(apart - sum) > clear) || !(std::abs(apart - gap) > clear))
      {
        if (!(apart > clear) || !Cleared(rounds, {.Length = clear, .Slack = Leeway(rounds)}, false))
        {
          return {};
        }
        const Vector toward = (second.Center - first.Center) / apart;
        const bool   inside = std::abs(apart - (first.Radius + second.Radius)) > clear;
        const Point  point  = inside && second.Radius > first.Radius
                                ? second.Center - (toward * second.Radius)
                                : first.Center + (toward * first.Radius);
        if (!(Swept(one, point) && Swept(other, point)))
        {
          return {.Kind = Tidings::Nothing};
        }
        Band band    = Rim(rounds, length);
        band.Parting = true;
        return band;
      }
      if (crossed)
      {
        const double cosine =
            ((inner * inner) + (outer * outer) - (apart * apart)) / (2.0 * inner * outer);
        if ((1.0 - (cosine * cosine)) > Met::Steep * Met::Steep)
        {
          return steep ? Riven(rounds, length) : Band{};
        }
        return Rim(rounds, length);
      }
      return {.Kind = Tidings::Nothing,
              .Lost = static_cast<int>(
                  !Cleared(rounds, {.Length = clear, .Slack = Leeway(rounds)}, true))};
    }

    /**
     *  What a band tells made marks: nothing, one touch at its point, or its crossings made exact,
     *  placed and put in order as Settled tells them.
     *
     *  @param  curves The two ready curves.
     *  @param  band   What the band tells.
     *  @param  ends   The carriers' two crossings, where the band has crossings.
     *  @return        The marks; nothing where Intersect must decide.
     */
    static auto Told(const Holds &curves, const Band &band,
                     const std::array<Chord::Settle, 2> &ends) -> std::optional<Abutment>
    {
      if (band.Kind == Tidings::Touch)
      {
        Abutment touch;
        touch.Marks[0] = {.At = band.At, .Kind = Contact::Touch};
        touch.Count    = 1;
        return touch;
      }
      if (band.Kind == Tidings::Nothing || band.Count == 0)
      {
        return Abutment{};
      }
      return Settled(curves, band, ends);
    }

    /**
     *  Whether Intersect's quick roots find a ready line's carrier and a ready arc's circle
     *  certainly apart, told from the discriminant worked out plainly: from the rounded way and
     *  gap, b^2 - a c is within 32 u a (G + r^2) of the true one, G = |P - O|^2, the roundings of
     *  each difference, product and sum taken; past that and Quadratic's 2^-94 a (7 G + 3 r^2) it
     *  lies past five times Intersect's bound.
     *
     *  @param  bow The line and the arc.
     *  @return     Whether it is told so; where not, Quadratic tells it.
     */
    static auto Parted(const Bow &bow) -> bool
    {
      const Arc   &arc    = Arced(*bow.Round);
      const Vector gap    = bow.Line->Ends[0] - arc.Center;
      const double b      = Dot(bow.Way, gap);
      const double away   = Dot(gap, gap);
      const double square = arc.Radius * arc.Radius;
      const double disc   = (b * b) - (bow.Square * (away - square));
      return std::abs(disc) > 0x1p-47 * bow.Square * (away + square);
    }

    /**
     *  Whether Intersect's quick roots find two ready arcs' circles certainly apart, told from big
     *  worked out plainly: from the rounded gap, 4 G r1^2 - k^2 is within 24 u S^2 of the true one,
     *  S = G + r1^2 + r2^2; past that and Circled's 2^-91 S^2 it lies past five times Intersect's
     *  bound.
     *
     *  @param  rounds The two ready arcs.
     *  @return        Whether it is told so; where not, Circled tells it.
     */
    static auto Parted(const Holds &rounds) -> bool
    {
      const Arc   &first  = Arced(rounds[0]);
      const Arc   &second = Arced(rounds[1]);
      const Vector gap    = second.Center - first.Center;
      const double apart  = Dot(gap, gap);
      const double inner  = first.Radius * first.Radius;
      const double outer  = second.Radius * second.Radius;
      const double chord  = (apart + inner) - outer;
      const double big    = (4.0 * apart * inner) - (chord * chord);
      const double size   = apart + inner + outer;
      return std::abs(big) > 0x1p-47 * size * size;
    }

    /**
     *  Crossings Roundly keeps, found as Finding finds them, made hits as Intersect's direct way
     *  makes them: at an end where an end kept lies within the tolerance of one, in order along the
     *  first curve, and where that is a whole circle, a hit at its seam first, as Seamed puts it,
     *  by the points Intersect places the hits at on it and its own point at its start.
     *
     *  @param  curves The two ready curves.
     *  @param  band   The ends kept, and the tolerance.
     *  @param  found  The crossings found, each on both curves.
     *  @return        The hits.
     */
    static auto Ruled(const Holds &curves, const Band &band, const Found &found) -> Abutment
    {
      Abutment told = found.Hits;
      for (Mark &mark : std::span{told.Marks}.first(told.Count))
      {
        int ended = 0;
        for (const Point end : std::span{band.Ends}.first(band.Ended))
        {
          ended = std::max(ended, static_cast<int>(Distance(end, mark.At) <= band.Length));
        }
        mark.Kind = ended > 0 ? Contact::End : Contact::Cross;
      }
      const Held &first  = curves[0];
      const int   circle = static_cast<int>(first.Round)
                         * static_cast<int>(Whole(first))
                         * static_cast<int>(told.Count == 2);
      if (circle == 0)
      {
        return told;
      }
      const std::array<const Curve *, 2> shapes = {first.Shape, curves[1].get().Shape};
      const Point                        seam   = At(*first.Shape, 0.0);
      std::array<int, 2>                 sealed = {};
      for (std::size_t index = 0; index < 2; ++index)
      {
        const Point placed = Crossings::Placed(shapes, {}, told.Marks.at(index).At)
                                 .value_or(Crossings::Crossed{})
                                 .On[0];
        sealed.at(index)   = static_cast<int>(Distance(placed, seam) <= band.Length);
      }
      const std::size_t move =
          static_cast<std::size_t>(sealed[1]) * (1U - static_cast<std::size_t>(sealed[0]));
      const std::array marks = told.Marks;
      told.Marks             = {marks.at(move), marks.at(1 - move)};
      return told;
    }

    /**
     *  A ready line and a ready arc Roundly gives up or finds nothing exact for, an end near the
     *  other curve, told as Shy tells them: every arc's end, a whole circle's seam too, farther
     *  than four tolerances from the line; each end of the line near the circle stopping short of
     *  it, the carrier's crossing nearer it past it, and a hit at that end where the arc's point
     *  facing it, as Facing places it, lies within the tolerance of the end, halfway between the
     *  two; and the crossings within the line and the sweep, found as Finding finds them. All in
     *  order along the first curve, by the places Intersect gives them where an arc is first.
     *
     *  @param  curves The ready line and arc.
     *  @param  bow    The line and the arc.
     *  @param  band   The carrier's crossings, and the tolerance.
     *  @return        The hits; nothing where Intersect must decide.
     */
    static auto Shied(const Holds &curves, const Bow &bow, const Band &band)
        -> std::optional<Abutment>
    {
      const Held           &line   = *bow.Line;
      const Held           &round  = *bow.Round;
      const Arc            &arc    = Arced(round);
      const bool            first  = &curves[0].get() == &line;
      const Ambit           span   = {.Length = 4.0 * band.Length, .Slack = bow.Slack};
      int                   lost   = static_cast<int>(!Cleared(bow, span, true));
      std::array<Mark, 4>   marks  = {};
      std::array<double, 4> places = {};
      std::size_t           count  = 0;
      for (std::size_t index = 0; index < 2; ++index)
      {
        const Point end = line.Ends.at(index);
        if (Beyond(end, round, span))
        {
          continue;
        }
        const std::size_t nearer =
            Distance(end, band.Carried[0]) <= Distance(end, band.Carried[1]) ? 0 : 1;
        const double along = band.Places.at(nearer);
        const auto   side  = static_cast<double>(index);
        lost += static_cast<int>(!(((1.0 - side) * -along) + (side * (along - 1.0)) > 0.0));
        const std::optional<double> facing =
            Earliest(Elliptic(arc), Heading(end - arc.Center).InRadians());
        const Point here  = At(*line.Shape, side);
        const Point there = At(*round.Shape, facing.value_or(0.0));
        if (facing && Distance(here, there) <= band.Length)
        {
          places.at(count) = first ? side : *facing;
          marks.at(count)  = {.At = Midpoint(here, there), .Kind = Contact::End};
          ++count;
        }
      }
      Band inner = {.Kind = Tidings::Crossings, .Length = band.Length};
      for (std::size_t side = 0; side < 2; ++side)
      {
        const double along = band.Places.at(side);
        const int    kept  = static_cast<int>(0.0 < along)
                           * static_cast<int>(along < 1.0)
                           * static_cast<int>(Swept(round, band.Carried.at(side)));
        if (kept > 0)
        {
          inner.Add(band.Carried.at(side));
        }
      }
      if (inner.Count > 0)
      {
        const Found found = Finding(curves, inner, Halved(bow, Quadratic(bow)));
        const std::array<const Curve *, 2> shapes  = {curves[0].get().Shape, curves[1].get().Shape};
        lost                                      += found.Lost + found.Off;
        for (std::size_t index = 0; index < found.Hits.Count; ++index)
        {
          const Point point = found.Hits.Marks.at(index).At;
          places.at(count) =
              first ? found.Places.at(index)
                    : Crossings::Placed(shapes, {}, point).value_or(Crossings::Crossed{}).Along[0];
          marks.at(count) = found.Hits.Marks.at(index);
          ++count;
        }
      }
      /**
       *  Every place apart from every other past Spread, so that each is in Intersect's order.
       */
      std::array<std::size_t, 4> order = {0, 1, 2, 3};
      std::ranges::sort(
          order, {}, [&places, count](std::size_t index) -> double
          { return index < count ? places.at(index) : std::numeric_limits<double>::infinity(); });
      for (std::size_t index = 1; index < count; ++index)
      {
        lost += static_cast<int>(
            !(places.at(order.at(index)) - places.at(order.at(index - 1)) > Spread));
      }
      if (lost + static_cast<int>(count > 2) > 0)
      {
        return std::nullopt;
      }
      Abutment sorted = {.Marks = {}, .Count = count};
      for (std::size_t index = 0; index < count; ++index)
      {
        sorted.Marks.at(index) = marks.at(order.at(index));
      }
      return sorted;
    }

    /**
     *  A ready line and a ready arc crossing at a sine past Steep, an end near the other curve, as
     *  Intersect tells them: Roundly's crossings, as Ruled makes them, where Intersect finds each
     *  exact one on both curves; else, an end near, as Shy tells them.
     *
     *  @param  curves The ready line and arc.
     *  @param  bow    The line and the arc.
     *  @param  band   What Roundly tells.
     *  @return        The hits; nothing where Intersect must decide.
     */
    static auto Pierce(const Holds &curves, const Bow &bow, const Band &band)
        -> std::optional<Abutment>
    {
      if (band.Kind == Tidings::Shy)
      {
        return Shied(curves, bow, band);
      }
      const Found found = Finding(curves, band, Halved(bow, Quadratic(bow)));
      if (found.Lost + (static_cast<int>(found.Off > 0) * static_cast<int>(!band.Near)) > 0)
      {
        return std::nullopt;
      }
      return found.Off == 0 ? std::optional<Abutment>{Ruled(curves, band, found)}
                            : Shied(curves, bow, band);
    }

    /**
     *  Two ready curves whose carriers Intersect asks to touch, where its quick roots may not find
     *  them apart: told as Intersect's Exactly tells them, by its own Crossing for a touch, at the
     *  places it asks at, its quick roots and its Tangency's short exact sums. Where they touch
     *  exactly, one touch, where Crossing places it on both curves; where they do not, Intersect
     *  goes on to its band, as Told tells it. The band is a touch or nothing: carriers not clearly
     *  apart lie within far less than a tolerance of touching.
     *
     *  @param  curves The two ready curves.
     *  @param  band   Intersect's band about where they come nearest, and where it asks for the
     *                 touch.
     *  @return        The marks; nothing where the band is not told.
     */
    static auto Kissed(const Holds &curves, const Band &band) -> std::optional<Abutment>
    {
      const std::array<const Curve *, 2> shapes = {curves[0].get().Shape, curves[1].get().Shape};
      const std::optional<Crossings::Crossed> touch = Crossings::Crossing(
          shapes, {.Along = band.Asked, .At = band.At, .Touch = true}, Crossings::Quickly(shapes));
      if (!touch)
      {
        return Told(curves, band, {});
      }
      Abutment touched;
      touched.Marks[0] = {.At = touch->At, .Kind = Contact::Touch};
      touched.Count    = 1;
      return touched;
    }

    /**
     *  A ready line and a ready arc, either first, told as Bent tells them: the carriers certainly
     *  apart where Intersect asks for their touch, plainly or else by Quadratic; their crossings
     *  from Quadratic where the band has crossings.
     *
     *  @param  curves The ready line and arc.
     *  @param  length The tolerance.
     *  @param  steep  Whether a steep crossing is told, or only bands.
     *  @return        The marks; nothing where Intersect must decide.
     */
    static auto Bowed(const Holds &curves, double length, bool steep) -> std::optional<Abutment>
    {
      const bool   first    = std::holds_alternative<Line>(*curves[0].get().Shape);
      const Held  &straight = curves.at(first ? 0 : 1);
      const Held  &round    = curves.at(first ? 1 : 0);
      const Arc   &arc      = Arced(round);
      const Vector way      = straight.Ends[1] - straight.Ends[0];
      const double size     = std::max({straight.Reach, std::abs(arc.Center.X) + arc.Radius,
                                        std::abs(arc.Center.Y) + arc.Radius});
      const Bow    bow      = {.Line   = &straight,
                               .Round  = &round,
                               .Way    = way,
                               .Square = Dot(way, way),
                               .Slack  = Geoxy::Leeway(size),
                               .Steep  = steep,
                               .First  = &curves[0].get()};
      const Band   band     = Bent(bow, length);
      if (static_cast<int>(band.Kind == Tidings::Unsure) + band.Lost > 0)
      {
        return std::nullopt;
      }
      if (band.Kind == Tidings::Steep || band.Kind == Tidings::Shy)
      {
        return Pierce(curves, bow, band);
      }
      std::optional<Quadric> shape;
      if (band.Parting && !Parted(bow))
      {
        shape = Quadratic(bow);
        if (!((std::abs(shape->Disc.High) * Lean::Tight) - shape->OffD
              > shape->Quick * Lean::Loose))
        {
          return Kissed(curves, band);
        }
      }
      const std::array<Chord::Settle, 2> ends = band.Count > 0
                                                  ? Halved(bow, shape.value_or(Quadratic(bow)))
                                                  : std::array<Chord::Settle, 2>{};
      return Told(curves, band, ends);
    }

    /**
     *  Two ready arcs told as Ringed tells them: the circles certainly apart where Intersect asks
     *  for their touch, plainly or else by Circled; their crossings from Circled where the band has
     *  crossings; and where they cross steeply, an end near, Circularly's crossings, as Ruled makes
     *  them, where Intersect finds each exact one on both arcs.
     *
     *  @param  rounds The two ready arcs.
     *  @param  length The tolerance.
     *  @param  steep  Whether a steep crossing is told, or only bands.
     *  @return        The marks; nothing where Intersect must decide.
     */
    static auto Rung(const Holds &rounds, double length, bool steep) -> std::optional<Abutment>
    {
      const Band band = Ringed(rounds, length, steep);
      if (static_cast<int>(band.Kind == Tidings::Unsure) + band.Lost > 0)
      {
        return std::nullopt;
      }
      if (band.Kind == Tidings::Steep)
      {
        const Found found = Finding(rounds, band, Halved(rounds, Circled(rounds)));
        return found.Lost + found.Off > 0 ? std::nullopt
                                          : std::optional<Abutment>{Ruled(rounds, band, found)};
      }
      std::optional<Circles> circle;
      if (band.Parting && !Parted(rounds))
      {
        circle = Circled(rounds);
        if (!((std::abs(circle->Big.High) * Lean::Tight) - circle->OffB
              > circle->Quick * Lean::Loose))
        {
          return Kissed(rounds, band);
        }
      }
      const std::array<Chord::Settle, 2> ends = band.Count > 0
                                                  ? Halved(rounds, circle.value_or(Circled(rounds)))
                                                  : std::array<Chord::Settle, 2>{};
      return Told(rounds, band, ends);
    }
  };

  auto Glancing::Hold(const Prepared &curve) -> Held
  {
    return {.Shape   = &curve.Shape(),
            .Ends    = curve.m_Ends,
            .Slip    = curve.m_Slip,
            .Sweep   = curve.m_Sweep,
            .Square  = curve.m_Square,
            .Under   = curve.m_Under,
            .Inverse = curve.m_Inverse,
            .Bounds  = curve.m_Box,
            .Reach   = curve.m_Reach,
            .Plain   = curve.m_Plain,
            .Round   = curve.m_Round};
  }

  auto Glancing::Hold(const Curve &curve, const std::array<Point, 2> &tips) -> Held
  {
    /**
     *  As Prepared makes them: its box Confines's, a line's between its ends, an arc's about its
     *  whole circle, and the largest magnitude of its edges; a line's way exactly in two doubles
     *  and its square, an arc's ends by the library's sine and cosine at its angles in radians and
     *  its radius squared, exactly in two doubles.
     */
    const std::optional<Box> box    = Confines(curve);
    const Box                bounds = box.value_or(Box{});
    const double             reach  = std::max({std::abs(bounds.Min.X), std::abs(bounds.Max.X),
                                                std::abs(bounds.Min.Y), std::abs(bounds.Max.Y)});
    const auto              *line   = std::get_if<Line>(&curve);
    Held                     held   = {.Shape  = &curve,
                                       .Bounds = bounds,
                                       .Reach  = reach,
                                       .Plain  = box.has_value() && std::isfinite(reach)};
    if (line != nullptr)
    {
      const Lean::Way    way    = Lean::Apart(line->End, line->Start);
      const Lean::Double square = Lean::Dot(way, way).Value;
      held.Slip                 = {.X = way[0].Low, .Y = way[1].Low};
      held.Square               = square.High;
      held.Under                = square.Low;
      held.Inverse              = square.High > 0.0 ? 1.0 / square.High : 0.0;
      held.Ends                 = {line->Start, line->End};
      return held;
    }
    if (!held.Plain)
    {
      return held;
    }
    const Arc   &round = std::get<Arc>(curve);
    const double sweep = round.Sweep.InRadians();
    const auto   cut   = [&round](Vector way) -> Point
    {
      return round.Center + (way * round.Radius);
    };
    const Twofold::Two square = Twofold::TwoProduct(round.Radius, round.Radius);
    held.Round  = round.Radius > 0.0 && std::abs(sweep) > 0.0 && std::abs(sweep) <= 2.0 * Pi;
    held.Square = square.High;
    held.Under  = square.Low;
    held.Sweep  = sweep;
    held.Ends   = {cut(Trig::Way(round.Start)), cut(Trig::Way(round.Start, round.Sweep))};
    held.Tips   = tips;
    held.Tipped = true;
    return held;
  }

  auto Glancing::Hinged(const Holds &curves, const Tolerance &sane) -> std::optional<Abutment>
  {
    const Held &one   = curves[0];
    const Held &other = curves[1];
    if (!one.Plain || !other.Plain)
    {
      return std::nullopt;
    }
    /**
     *  The length Intersect's Reached gives two plain curves, and nothing where their boxes lie
     *  farther than four of it apart, as it tells them; then Steeply, as Intersect asks it within
     *  Coarse and Fine.
     */
    const double length = std::max(sane.Length, 16.0 * Above(std::max(one.Reach, other.Reach)));
    const Box   &box    = one.Bounds;
    const Box   &next   = other.Bounds;
    const double clear  = 4.0 * length;
    if (box.Min.X - next.Max.X    > clear
        || next.Min.X - box.Max.X > clear
        || box.Min.Y - next.Max.Y > clear
        || next.Min.Y - box.Max.Y > clear)
    {
      return Abutment{};
    }
    const bool line   = std::holds_alternative<Line>(*one.Shape);
    const bool across = std::holds_alternative<Line>(*other.Shape);
    if (length > Coarse || length < Fine || (line && across))
    {
      return std::nullopt;
    }
    if (line || across)
    {
      return (line ? other : one).Round ? Bands::Bowed(curves, length, true) : std::nullopt;
    }
    return one.Round && other.Round ? Bands::Rung(curves, length, true) : std::nullopt;
  }

  auto Glancing::Met(const Prepared::Pair &curves, const Tolerance &usable)
      -> std::optional<Abutment>
  {
    const Prepared &one    = curves[0];
    const Prepared &other  = curves[1];
    const bool      line   = std::holds_alternative<Line>(one.Shape());
    const bool      across = std::holds_alternative<Line>(other.Shape());
    if (!one.m_Plain || !other.m_Plain || (line && across))
    {
      return std::nullopt;
    }
    const Held  first  = Hold(one);
    const Held  second = Hold(other);
    const Holds held   = {first, second};
    if (line || across)
    {
      return (line ? other : one).m_Round ? Bands::Bowed(held, usable.Length, false) : std::nullopt;
    }
    return one.m_Round && other.m_Round ? Bands::Rung(held, usable.Length, false) : std::nullopt;
  }

  auto Glancing::Steeply(const Prepared::Pair &curves, const Tolerance &usable)
      -> std::optional<Abutment>
  {
    const Prepared &one    = curves[0];
    const Prepared &other  = curves[1];
    const bool      line   = std::holds_alternative<Line>(one.Shape());
    const bool      across = std::holds_alternative<Line>(other.Shape());
    if (!one.m_Plain || !other.m_Plain || (line && across))
    {
      return std::nullopt;
    }
    const Held  first  = Hold(one);
    const Held  second = Hold(other);
    const Holds held   = {first, second};
    if (line || across)
    {
      return (line ? other : one).m_Round ? Bands::Bowed(held, usable.Length, true) : std::nullopt;
    }
    return one.m_Round && other.m_Round ? Bands::Rung(held, usable.Length, true) : std::nullopt;
  }
} // namespace Geoxy
