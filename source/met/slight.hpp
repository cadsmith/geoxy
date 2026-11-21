/***************************************************************************************************
 *  @file      source/met/slight.hpp
 *  @brief     A polyline's remote arcs, too slight to be worked out through their centers, met
 *             through their chords and told again on the arcs themselves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_SLIGHT_HPP
#define GEOXY_MET_SLIGHT_HPP

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <optional>

namespace Geoxy::Met
{
  /**
   *  Whether a polyline has a remote arc, as Remote tells one.
   *
   *  @param  polyline The polyline.
   *  @return          Whether it has one.
   */
  auto Slightish(const Polyline &polyline) -> bool;

  /**
   *  Whether a curve is a polyline with a remote arc.
   *
   *  @param  curve The curve.
   *  @return       Whether it is.
   */
  auto Slightish(const Curve &curve) -> bool;

  /**
   *  Whether either of two curves has a remote arc, where that is not known already: while a
   *  Vetting stands, the curves being met were looked at and have none.
   *
   *  @param  first  One curve.
   *  @param  second The other.
   *  @return        Whether either has one; never while a Vetting stands.
   */
  auto Remotely(const Curve &first, const Curve &second) -> bool;

  /**
   *  The curves being met told to have no remote arc, for as long as this stands; not while the
   *  thorough way is meeting curves that have some, which it was asked to by Slightly.
   */
  class Vetting
  {
  public:
    Vetting();
    ~Vetting();
    Vetting(const Vetting &)                     = delete;
    auto operator=(const Vetting &) -> Vetting & = delete;
    Vetting(Vetting &&)                          = delete;
    auto operator=(Vetting &&) -> Vetting &      = delete;

  private:
    bool m_Before = false; /**< Whether one stood already. */
  };

  /**
   *  A polyline with every remote arc made its chord, to stand for it.
   */
  struct Standing
  {
    Curve  Shape   = Line{}; /**< The polyline with its chords.                              */
    double Deepest = 0.0;    /**< How far the remote arc farthest from its chord stands off. */
  };

  /**
   *  What stands for a curve with remote arcs.
   *
   *  @param  curve The curve.
   *  @return       The polyline with its chords; nothing for a curve with no remote arc.
   */
  auto Stood(const Curve &curve) -> std::optional<Standing>;

  /**
   *  How two curves with remote arcs are met.
   */
  struct Chording
  {
    Tolerance Sane   = {};  /**< The tolerance, as Sane makes it. */
    double    Length = 0.0; /**< The length used: the tolerance's, no finer than the curves'
                                 grain.                                                      */
    void (*Intersect)(const void *, const Curves &, const Tolerance &, Intersections &) =
        nullptr; /**< What meets the curves that stand for them, given the context. */
    void (*Deep)(const Curves &, const Tolerance &, Intersections &) =
        nullptr; /**< What meets curves with remote arcs the thorough way.          */
    const void *Context                  = nullptr; /**< What Intersect is given first.          */
    std::array<const Standing *, 2> Kept = {nullptr, nullptr}; /**< What stands for each
                                                                    curve, where it is kept. */
  };

  /**
   *  Two curves met where one has remote arcs.
   *
   *  Every remote arc's chord stands for it, which every direct way takes. An arc within the
   *  tolerance of its chord is, within the tolerance, the chord. Where one stands off by more, the
   *  curves are met again with the tolerance widened by as much: where that finds, on such arcs,
   *  nothing but the same plain crossings, each steep beside the arc's own turn and clear of its
   *  ends, each is one crossing of the arc itself and nothing else of it is near. A crossing found
   *  so is put where the curves themselves cross, and a hit's point on a chord is the arc's own at
   *  its place. Where that is not found, the thorough way meets the curves, each such arc a piece.
   *
   *  @param  curves The two curves.
   *  @param  asked  The tolerance, the length used and what meets the curves.
   *  @param  into   Where what they meet in goes.
   *  @return        Whether they were met so: whether either has a remote arc.
   */
  auto Slightly(const Curves &curves, const Chording &asked, Intersections &into) -> bool;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_SLIGHT_HPP */
