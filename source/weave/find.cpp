/***************************************************************************************************
 *  @file      source/weave/find.cpp
 *  @brief     Every crossing and shared stretch of an arrangement's curves, found from the pairs
 *             that come near.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "weave/find.hpp"

#include "abut.hpp"
#include "crossing.hpp"
#include "graze.hpp"
#include "sane.hpp"
#include "scratch.hpp"
#include "weave/cut.hpp"
#include "weave/found.hpp"
#include "weave/pairs.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  namespace
  {
    /**
     *  How many tolerances two curves' boxes are grown by before they are asked whether they come
     *  near: four times the sixteen within which Intersect tells two curves apart from their boxes.
     */
    constexpr double Reaching = 64.0;

    /**
     *  A held line's or arc's ends: a line's own, an arc's as At gives them.
     *
     *  @param  held The curve, held.
     *  @return      Its start and end.
     */
    auto Tipped(const Glancing::Held &held) -> const std::array<Point, 2> &
    {
      return held.Tipped ? held.Tips : held.Ends;
    }

    /**
     *  Whether two held curves are joined end to end, as a polyline's segments are: an end of one
     *  within a length of an end of the other, by their largest coordinate gap.
     *
     *  @param  holds  The two curves, held.
     *  @param  length The length.
     *  @return        Whether they are.
     */
    auto Kinked(const std::array<const Glancing::Held *, 2> &holds, double length) -> bool
    {
      const std::array<Point, 2> &one   = Tipped(*holds[0]);
      const std::array<Point, 2> &other = Tipped(*holds[1]);
      return std::ranges::any_of(
          one, [&other, length](Point end) -> bool
          { return Largest(end - other[0]) <= length || Largest(end - other[1]) <= length; });
    }

    /**
     *  The pairs of curves met, each curve made ready, or held for the bands, the first time a pair
     *  asks for either: on a polygon most sides meet only their neighbours, at joints told from the
     *  lines alone, and a bulged one's arcs theirs, told by the bands.
     */
    class Meetings
    {
    public:
      /**
       *  Curves of which none is yet made ready.
       *
       *  @param sources   The lines and arcs.
       *  @param holds     Each arc held, or none for no arc; the rest are held as asked.
       *  @param tolerance The tolerance.
       */
      Meetings(std::span<const Curve> sources, List<Glancing::Held> holds,
               const Tolerance &tolerance)
          : m_Sources(sources), m_Holds(std::move(holds)), m_Made(sources.size(), sources.size()),
            m_Tolerance(tolerance), m_Sane(Sane(tolerance))
      {
        m_Readies.reserve(sources.size());
      }

      /**
       *  What two curves meet in, into a marking, each end hit at where the curves cross:
       *  Intersect's end hit stands at the end lying within the tolerance of the other curve,
       *  however far along the band about a shallow crossing the two cross, and the arrangement's
       *  vertex is where they cross, nearest it on both, where they do.
       *
       *  @param pair    The two curves.
       *  @param kind    Their pair's kind.
       *  @param marking The marking, whatever it held replaced.
       */
      void Met(const Couple &pair, std::size_t kind, Marking &marking)
      {
        Gathered(pair, kind, marking);
        const std::array<const Curve *, 2> curves = {&m_Sources[pair[0]], &m_Sources[pair[1]]};
        for (Mark &mark : marking.Marks)
        {
          if (mark.Kind != Contact::End)
          {
            continue;
          }
          const auto place = [&mark](const Curve &curve) -> double
          {
            const std::vector<Place> feet = Perpendiculars(curve, mark.At);
            return feet.empty() ? 0.0 : feet.front().Along;
          };
          const std::optional<Crossings::Crossed> crossing =
              Crossings::Crossing(curves, {place(*curves[0]), place(*curves[1])}, mark.At);
          if (crossing)
          {
            mark.At = crossing->At;
          }
        }
      }

    private:
      /**
       *  What two curves meet in, into a marking: a line and an arc, or two arcs, as Hinged tells
       *  them, which answers as Intersect does where it answers, neither made ready; the rest as
       *  Intersect tells them made ready.
       *
       *  @param pair    The two curves.
       *  @param kind    Their pair's kind.
       *  @param marking The marking, whatever it held replaced.
       */
      void Gathered(const Couple &pair, std::size_t kind, Marking &marking)
      {
        const bool                    banded = kind >= Kinds - 2;
        const std::optional<Abutment> told =
            banded ? Glancing::Hinged({Holding(pair[0]), Holding(pair[1])}, m_Sane) : std::nullopt;
        /**
         *  A joint the bands leave is one Intersect works out the thorough way, its passes before
         *  it telling nothing: met as curves, as Intersect made ready meets them, with none made
         *  ready.
         */
        if (!told && banded && Kinked({&m_Holds[pair[0]], &m_Holds[pair[1]]}, m_Sane.Length))
        {
          Geoxy::Intersect(Curves{m_Sources[pair[0]], m_Sources[pair[1]]}, m_Tolerance, m_Meeting);
          marking.Marks.clear();
          for (const Hit &hit : m_Meeting.Hits)
          {
            marking.Marks.push_back({.At = hit.At, .Kind = hit.Kind});
          }
          marking.Overlaps.assign(m_Meeting.Overlaps.begin(), m_Meeting.Overlaps.end());
          marking.Spent = m_Meeting.Spent;
          return;
        }
        if (!told)
        {
          Intersect(Prepared::Pair{Readied(pair[0]), Readied(pair[1])}, m_Tolerance, marking);
          return;
        }
        const std::span<const Mark>    marks  = std::span{told->Marks}.first(told->Count);
        const std::span<const Overlap> shared = std::span{&told->Stretch, told->Stretches};
        marking.Marks.assign(marks.begin(), marks.end());
        marking.Overlaps.assign(shared.begin(), shared.end());
        marking.Spent = false;
      }

      /**
       *  A curve held for the bands, held the first time it is asked for.
       *
       *  @param  index The curve's index.
       *  @return       It, held.
       */
      auto Holding(std::size_t index) -> const Glancing::Held &
      {
        if (m_Holds[index].Shape == nullptr)
        {
          m_Holds[index] = Glancing::Hold(m_Sources[index], {});
        }
        return m_Holds[index];
      }

      /**
       *  A curve made ready, made the first time it is asked for.
       *
       *  @param  index The curve's index.
       *  @return       It, ready.
       */
      auto Readied(std::size_t index) -> const Prepared &
      {
        if (m_Made[index] == m_Sources.size())
        {
          m_Made[index] = m_Readies.size();
          m_Readies.emplace_back(m_Sources[index]);
        }
        return m_Readies[m_Made[index]];
      }

      std::span<const Curve> m_Sources; /**< The lines and arcs.                   */
      List<Glancing::Held>   m_Holds;   /**< Each curve held, once it is.          */
      List<Prepared>         m_Readies; /**< The curves made ready, in that order. */
      List<std::size_t>      m_Made;    /**< Each curve's place among them, or past
                                             the last curve where it is not made. */
      Tolerance     m_Tolerance;        /**< The tolerance.                         */
      Tolerance     m_Sane;             /**< The tolerance, as Sane makes it.       */
      Intersections m_Meeting = {};     /**< What two curves met as curves meet in. */
    };

    /**
     *  What Find knows of each curve before it lays its vertices: its stance, whether it is busy,
     *  and its ends as At gives them.
     */
    struct Facts
    {
      List<std::uint8_t> Stances = {};      /**< Each curve's stance. */
      List<std::uint8_t> Busy    = {};      /**< Whether each is in a pair found,
                                                 or in a joint told with a curve
                                                 not next to it.                 */
      List<std::array<Point, 2>> Ends = {}; /**< Each curve's ends. */
    };

    /**
     *  Which joints of two lines in a row are lone: each line meeting no curve but its neighbours,
     *  at their joints, so neither busy, the first's end the very point of the next's start, and
     *  both lines longer than the length within which vertices merge.
     *
     *  Such a joint, made a vertex, would stand for itself among the vertices merged within that
     *  length, and merge with its twin alone; no edge but the two lines' would reach it. Every
     *  vertex lies within the length of a point of a curve it is of: an end, a mark of a pair met,
     *  or a shared stretch's end. One within the length of the joint, of a curve neither line, has
     *  its curve's box, grown by Reaching times as much, meet the box of each line, grown alike, so
     *  the curve is in a pair with each, as Woven finds every such pair; the curve is next to one
     *  of them at most, and in a pair with the other it makes that one busy. The lines' own other
     *  ends lie farther than the length, and they have no marks.
     *
     *  @param  facts  The curves' stances, whether each is busy, and their ends.
     *  @param  length The length within which vertices merge.
     *  @return        Whether each curve's end and the next's start are a lone joint.
     */
    auto Lonely(const Facts &facts, double length) -> List<std::uint8_t>
    {
      const std::size_t  count = facts.Stances.size();
      List<std::uint8_t> strung(count, 0);
      for (std::size_t curve = 0; curve + 1 < count; ++curve)
      {
        const std::array<Point, 2> &one   = facts.Ends[curve];
        const std::array<Point, 2> &other = facts.Ends[curve + 1];
        if (facts.Stances[curve] < 3
            && facts.Stances[curve + 1] < 3
            && facts.Busy[curve] == 0
            && facts.Busy[curve + 1] == 0
            && Twinned(one[1], other[0]))
        {
          const double shorter = std::min(Distance(one[0], one[1]), Distance(other[0], other[1]));
          strung[curve]        = static_cast<std::uint8_t>(shorter > length);
        }
      }
      return strung;
    }

    /**
     *  What tells the joints of lines in a row: each line's limb, whether Seamless tells each
     *  line's joint with the next, and the tolerance as Sane makes it.
     */
    struct Seams
    {
      std::span<const Limb>         Limbs  = {}; /**< Each line's limb.           */
      std::span<const std::uint8_t> Joined = {}; /**< Whether each joint is told. */
      Tolerance                     Sane   = {}; /**< The tolerance, made sane.   */
    };

    /**
     *  The pairs of curves whose boxes overlap that are not joints told, each curve in one, or in a
     *  joint told with any curve but the next or the last, made busy: so that a curve in none meets
     *  nothing but its neighbours, at their joints.
     *
     *  @param  boxes The curves' boxes, two or more.
     *  @param  facts The curves' stances, and whether each is busy, marked.
     *  @param  seams What tells the joints.
     *  @return       The pairs, in no order.
     */
    auto Woof(std::span<const Box> boxes, Facts &facts, const Seams &seams) -> List<Couple>
    {
      List<Couple> found;
      found.reserve(boxes.size() / 8);
      Woven(boxes,
            [&facts, &seams, &found](std::size_t first, std::size_t second) -> void
            {
              if (second == first + 1 && seams.Joined[first] != 0)
              {
                return;
              }
              const bool lines =
                  Couplings.at(facts.Stances[first]).at(facts.Stances[second]) < Kinds - 2;
              const bool told =
                  lines && Abutting::Jointed({seams.Limbs[first], seams.Limbs[second]}, seams.Sane);
              if (!told || second != first + 1)
              {
                facts.Busy[first]  = 1;
                facts.Busy[second] = 1;
              }
              if (!told)
              {
                found.push_back({first, second});
              }
            });
      return found;
    }
  } // namespace

  auto Find(std::span<const Curve> sources, const Tolerance &tolerance, const Compass &compass)
      -> Findings
  {
    const std::size_t    count = sources.size();
    const double         pad   = Reaching * compass.Reach;
    const Tolerance      sane  = Sane(tolerance);
    Cutter               cutter{count};
    List<Box>            boxes(count);
    Facts                facts = {.Stances = List<std::uint8_t>(count),
                                  .Busy    = List<std::uint8_t>(count, 0),
                                  .Ends    = List<std::array<Point, 2>>(count)};
    List<Limb>           limbs(count);
    List<std::uint8_t>   joined(count, 0);
    List<Glancing::Held> holds;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Curve &source  = sources[index];
      const auto  *line    = std::get_if<Line>(&source);
      facts.Stances[index] = Stance(source);
      if (line == nullptr)
      {
        const Arc                 &arc  = std::get<Arc>(source);
        const std::array<Point, 2> ends = {At(arc, 0.0), At(arc, 1.0)};
        facts.Ends[index]               = ends;
        boxes[index]                    = Padded(Snug(arc, ends), pad);
        holds.resize(count);
        holds[index] = Glancing::Hold(source, ends);
        continue;
      }
      /**
       *  A line's box is Snug's: its ends as At gives them are its own but for the sign of a zero,
       *  which a box padded by a length keeps nothing of.
       */
      const Point start = line->Start;
      const Point end   = line->End;
      facts.Ends[index] = {Lerp(start, end, 0.0), Lerp(start, end, 1.0)};
      boxes[index] = Padded({.Min = {.X = std::min(start.X, end.X), .Y = std::min(start.Y, end.Y)},
                             .Max = {.X = std::max(start.X, end.X), .Y = std::max(start.Y, end.Y)}},
                            pad);
      limbs[index] = Abutting::Limbed(*line);
      if (index > 0 && limbs[index - 1].Shape != nullptr)
      {
        joined[index - 1] =
            static_cast<std::uint8_t>(Seamless(limbs[index - 1], limbs[index], sane));
      }
    }
    const List<Couple> found =
        count > 1 ? Woof(boxes, facts, {.Limbs = limbs, .Joined = joined, .Sane = sane})
                  : List<Couple>{};
    /**
     *  The vertices laid in the curves' order, but a lone joint's where chained.
     */
    List<std::uint8_t> strung =
        compass.Chained ? Lonely(facts, compass.Length) : List<std::uint8_t>(count, 0);
    for (std::size_t index = 0; index < count; ++index)
    {
      cutter.Placed(index, facts.Ends[index],
                    {index > 0 && strung[index - 1] != 0, strung[index] != 0});
      if (const auto *arc = std::get_if<Arc>(&sources[index]))
      {
        cutter.Halved(index, *arc);
      }
    }
    const std::array<List<Couple>, Kinds> kinds = Coupled(found, facts.Stances);
    Meetings                              meetings{sources, std::move(holds), tolerance};
    Marking                               marking;
    for (std::size_t kind = 0; kind < Kinds; ++kind)
    {
      for (const Couple &pair : kinds.at(kind))
      {
        meetings.Met(pair, kind, marking);
        cutter.Record(pair, marking, sources);
      }
    }
    Findings findings = cutter.Found();
    if (compass.Chained)
    {
      findings.Strung = std::move(strung);
    }
    return findings;
  }
} // namespace Geoxy::Weave
