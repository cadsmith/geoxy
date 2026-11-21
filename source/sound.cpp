/***************************************************************************************************
 *  @file      source/sound.cpp
 *  @brief     Whether every number a curve is held by is finite, each number's product with zero
 *             summed, so that one test tells of them all.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "sound.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <cmath>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  A sum of numbers each times zero: zero while every one is finite, not a number once one is
     *  not, as an infinity or not a number times zero is not a number.
     */
    class Naughts
    {
    public:
      /**
       *  A number added.
       *
       *  @param value The number.
       */
      void Add(double value) { m_Sum += value * 0.0; }

      /**
       *  A point's coordinates added.
       *
       *  @param point The point.
       */
      void Add(Point point)
      {
        Add(point.X);
        Add(point.Y);
      }

      /**
       *  A vector's components added.
       *
       *  @param vector The vector.
       */
      void Add(Vector vector)
      {
        Add(vector.X);
        Add(vector.Y);
      }

      /**
       *  An angle added, in the unit it is held in or in degrees.
       *
       *  @param angle The angle.
       */
      void Add(Angle angle) { Add(angle.InDegrees()); }

      /**
       *  A curve's numbers added.
       *
       *  @param line The curve.
       */
      void Add(const Line &line)
      {
        Add(line.Start);
        Add(line.End);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param ray The curve.
       */
      void Add(const Ray &ray)
      {
        Add(ray.Start);
        Add(ray.Way);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param xline The curve.
       */
      void Add(const Xline &xline)
      {
        Add(xline.Through);
        Add(xline.Way);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param arc The curve.
       */
      void Add(const Arc &arc)
      {
        Add(arc.Center);
        Add(arc.Radius);
        Add(arc.Start);
        Add(arc.Sweep);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param ellipse The curve.
       */
      void Add(const Ellipse &ellipse)
      {
        Add(ellipse.Map.X);
        Add(ellipse.Map.Y);
        Add(ellipse.Map.Origin);
        Add(ellipse.Start);
        Add(ellipse.Sweep);
      }

      /**
       *  A curve's numbers added, its vertices' coordinates and bulges each in a sum of its own, so
       *  that no sum waits on another: each is 0 or not a number as the whole one would be.
       *
       *  @param polyline The curve.
       */
      void Add(const Polyline &polyline)
      {
        double across = 0.0;
        double along  = 0.0;
        double bulges = 0.0;
        for (const Vertex &vertex : polyline.Vertices)
        {
          across += vertex.Position.X * 0.0;
          along  += vertex.Position.Y * 0.0;
          bulges += vertex.Bulge * 0.0;
        }
        Add(across);
        Add(along);
        Add(bulges);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param spline The curve.
       */
      void Add(const Spline &spline)
      {
        for (const double knot : spline.Knots)
        {
          Add(knot);
        }
        for (const Point pole : spline.Poles)
        {
          Add(pole);
        }
        for (const double weight : spline.Weights)
        {
          Add(weight);
        }
      }

      /**
       *  A curve's numbers added.
       *
       *  @param quadratic The curve.
       */
      void Add(const Bezier::Quadratic &quadratic)
      {
        Add(quadratic.Start);
        Add(quadratic.Control);
        Add(quadratic.End);
      }

      /**
       *  A curve's numbers added.
       *
       *  @param cubic The curve.
       */
      void Add(const Bezier::Cubic &cubic)
      {
        Add(cubic.Start);
        Add(cubic.Out);
        Add(cubic.In);
        Add(cubic.End);
      }

      /**
       *  Whether every number added was finite.
       *
       *  @return Whether it was.
       */
      [[nodiscard]] auto Finite() const -> bool { return !std::isnan(m_Sum); }

    private:
      double m_Sum = 0.0; /**< The sum. */
    };
  } // namespace

  auto Sound(const Curve &curve) -> bool
  {
    Naughts naughts;
    std::visit([&naughts](const auto &shape) -> void { naughts.Add(shape); }, curve);
    return naughts.Finite();
  }
} // namespace Geoxy
