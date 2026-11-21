/***************************************************************************************************
 *  @file      tests/support/text.cpp
 *  @brief     Curves written as a line of text, and read back.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "support/text.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <charconv>
#include <cstddef>
#include <format>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The words of a text written one after another, a space between each.
     */
    class Writer
    {
    public:
      /**
       *  Writes a word.
       *
       *  @param word The word.
       */
      void Word(std::string_view word)
      {
        if (!m_Text.empty())
        {
          m_Text += ' ';
        }
        m_Text += word;
      }

      /**
       *  Writes a number, the shortest decimal that reads back as it.
       *
       *  @param number The number.
       */
      void Number(double number) { Word(std::format("{}", number)); }

      /**
       *  Writes an angle in its unit: radians as a number, degrees as a number and a d.
       *
       *  @param angle The angle.
       */
      void Turn(Angle angle)
      {
        if (angle.HeldInDegrees())
        {
          Word(std::format("{}d", angle.InDegrees()));
          return;
        }
        Number(angle.InRadians());
      }

      /**
       *  Writes a count.
       *
       *  @param count The count.
       */
      void Count(std::size_t count) { Word(std::format("{}", count)); }

      /**
       *  Writes a point, its x then its y.
       *
       *  @param point The point.
       */
      void At(Point point)
      {
        Number(point.X);
        Number(point.Y);
      }

      /**
       *  Writes a displacement, its x then its y.
       *
       *  @param vector The displacement.
       */
      void Way(Vector vector)
      {
        Number(vector.X);
        Number(vector.Y);
      }

      /**
       *  The text written.
       *
       *  @return The text.
       */
      [[nodiscard]] auto Text() -> std::string { return std::move(m_Text); }

    private:
      std::string m_Text; /**< The text so far. */
    };

    /**
     *  The words of a text read one after another; the first that fails marks the text not a curve.
     */
    class Reader
    {
    public:
      /**
       *  Reads a text.
       *
       *  @param text The text.
       */
      explicit Reader(std::string_view text) : m_Text(text) {}

      /**
       *  The next word.
       *
       *  @return The word, or none, which fails the text, past the last.
       */
      auto Word() -> std::string_view
      {
        Skip();
        const std::size_t      end  = m_Text.find_first_of(Blanks);
        const std::string_view word = m_Text.substr(0, end);
        m_Text.remove_prefix(word.size());
        m_Failed = m_Failed || word.empty();
        return word;
      }

      /**
       *  The next word as a number.
       *
       *  @return The number, or 0, which fails the text, when it is not one.
       */
      auto Number() -> double
      {
        double number = 0.0;
        Parse(Word(), number);
        return number;
      }

      /**
       *  The next word as an angle: radians, or degrees when it ends in d.
       *
       *  @return The angle, or none, which fails the text, when it is not one.
       */
      auto Turn() -> Angle
      {
        std::string_view word   = Word();
        double           number = 0.0;
        if (word.ends_with('d'))
        {
          word.remove_suffix(1);
          Parse(word, number);
          return Degrees(number);
        }
        Parse(word, number);
        return Radians(number);
      }

      /**
       *  The next word as a count of what follows, no more than the words left.
       *
       *  @return The count, or 0, which fails the text, when it is not one or runs past the words.
       */
      auto Count() -> std::size_t
      {
        std::size_t count = 0;
        Parse(Word(), count);
        if (count > Left())
        {
          m_Failed = true;
          return 0;
        }
        return count;
      }

      /**
       *  The next word as a whole number.
       *
       *  @return The number, or 0, which fails the text, when it is not one.
       */
      auto Whole() -> int
      {
        int whole = 0;
        Parse(Word(), whole);
        return whole;
      }

      /**
       *  The next two words as a point.
       *
       *  @return The point.
       */
      auto At() -> Point { return {.X = Number(), .Y = Number()}; }

      /**
       *  The next two words as a displacement.
       *
       *  @return The displacement.
       */
      auto Way() -> Vector { return {.X = Number(), .Y = Number()}; }

      /**
       *  Fails the text.
       */
      void Fail() { m_Failed = true; }

      /**
       *  Whether every word read was what it had to be, and none is left.
       *
       *  @return Whether the text is a curve.
       */
      auto Done() -> bool
      {
        Skip();
        return !m_Failed && m_Text.empty();
      }

    private:
      static constexpr std::string_view Blanks = " \t\r\n"; /**< What lies between words. */

      /**
       *  Parses a word wholly, failing the text when it is not the kind asked.
       *
       *  @tparam ValueT What the word is.
       *  @param  word   The word.
       *  @param  value  Where it goes.
       */
      template <typename ValueT>
      void Parse(std::string_view word, ValueT &value)
      {
        const auto [end, error] =
            std::from_chars(std::to_address(word.begin()), std::to_address(word.end()), value);
        m_Failed = m_Failed || error != std::errc{} || end != std::to_address(word.end());
      }

      /**
       *  Steps past the blanks before the next word.
       */
      void Skip()
      {
        const std::size_t start = m_Text.find_first_not_of(Blanks);
        m_Text.remove_prefix(start == std::string_view::npos ? m_Text.size() : start);
      }

      /**
       *  The number of words left.
       *
       *  @return The count.
       */
      auto Left() -> std::size_t
      {
        Reader      copy  = *this;
        std::size_t count = 0;
        while (!copy.Word().empty())
        {
          ++count;
        }
        return count;
      }

      std::string_view m_Text;           /**< What is left to read.                     */
      bool             m_Failed = false; /**< Whether a word was not what it had to be. */
    };

    /**
     *  Writes a line.
     *
     *  @param out  The text.
     *  @param line The line.
     */
    void Write(Writer &out, const Line &line)
    {
      out.Word("line");
      out.At(line.Start);
      out.At(line.End);
    }

    /**
     *  Writes a ray.
     *
     *  @param out The text.
     *  @param ray The ray.
     */
    void Write(Writer &out, const Ray &ray)
    {
      out.Word("ray");
      out.At(ray.Start);
      out.Way(ray.Way);
    }

    /**
     *  Writes an xline.
     *
     *  @param out   The text.
     *  @param xline The xline.
     */
    void Write(Writer &out, const Xline &xline)
    {
      out.Word("xline");
      out.At(xline.Through);
      out.Way(xline.Way);
    }

    /**
     *  Writes an arc.
     *
     *  @param out The text.
     *  @param arc The arc.
     */
    void Write(Writer &out, const Arc &arc)
    {
      out.Word("arc");
      out.At(arc.Center);
      out.Number(arc.Radius);
      out.Turn(arc.Start);
      out.Turn(arc.Sweep);
    }

    /**
     *  Writes an ellipse.
     *
     *  @param out     The text.
     *  @param ellipse The ellipse.
     */
    void Write(Writer &out, const Ellipse &ellipse)
    {
      out.Word("ellipse");
      out.Way(ellipse.Map.X);
      out.Way(ellipse.Map.Y);
      out.At(ellipse.Map.Origin);
      out.Turn(ellipse.Start);
      out.Turn(ellipse.Sweep);
    }

    /**
     *  Writes a polyline.
     *
     *  @param out      The text.
     *  @param polyline The polyline.
     */
    void Write(Writer &out, const Polyline &polyline)
    {
      out.Word("polyline");
      out.Word(polyline.Closed ? "closed" : "open");
      out.Count(polyline.Vertices.size());
      for (const Vertex &vertex : polyline.Vertices)
      {
        out.At(vertex.Position);
        out.Number(vertex.Bulge);
      }
    }

    /**
     *  Writes a spline.
     *
     *  @param out    The text.
     *  @param spline The spline.
     */
    void Write(Writer &out, const Spline &spline)
    {
      out.Word("spline");
      out.Word(std::format("{}", spline.Degree));
      out.Count(spline.Knots.size());
      for (const double knot : spline.Knots)
      {
        out.Number(knot);
      }
      out.Count(spline.Poles.size());
      for (const Point pole : spline.Poles)
      {
        out.At(pole);
      }
      out.Count(spline.Weights.size());
      for (const double weight : spline.Weights)
      {
        out.Number(weight);
      }
    }

    /**
     *  Writes a quadratic Bézier curve.
     *
     *  @param out       The text.
     *  @param quadratic The curve.
     */
    void Write(Writer &out, const Bezier::Quadratic &quadratic)
    {
      out.Word("quadratic");
      out.At(quadratic.Start);
      out.At(quadratic.Control);
      out.At(quadratic.End);
    }

    /**
     *  Writes a cubic Bézier curve.
     *
     *  @param out   The text.
     *  @param cubic The curve.
     */
    void Write(Writer &out, const Bezier::Cubic &cubic)
    {
      out.Word("cubic");
      out.At(cubic.Start);
      out.At(cubic.Out);
      out.At(cubic.In);
      out.At(cubic.End);
    }

    /**
     *  Reads a polyline's words after its kind.
     *
     *  @param  words The text.
     *  @return       The polyline.
     */
    auto ReadPolyline(Reader &words) -> Polyline
    {
      const std::string_view closure = words.Word();
      if (closure != "open" && closure != "closed")
      {
        words.Fail();
      }
      Polyline          polyline = {.Vertices = {}, .Closed = closure == "closed"};
      const std::size_t count    = words.Count();
      for (std::size_t index = 0; index < count; ++index)
      {
        const Point position = words.At();
        polyline.Vertices.push_back({.Position = position, .Bulge = words.Number()});
      }
      return polyline;
    }

    /**
     *  Reads a spline's words after its kind.
     *
     *  @param  words The text.
     *  @return       The spline.
     */
    auto ReadSpline(Reader &words) -> Spline
    {
      Spline            spline = {.Degree = words.Whole(), .Knots = {}, .Poles = {}, .Weights = {}};
      const std::size_t knots  = words.Count();
      for (std::size_t index = 0; index < knots; ++index)
      {
        spline.Knots.push_back(words.Number());
      }
      const std::size_t poles = words.Count();
      for (std::size_t index = 0; index < poles; ++index)
      {
        spline.Poles.push_back(words.At());
      }
      const std::size_t weights = words.Count();
      for (std::size_t index = 0; index < weights; ++index)
      {
        spline.Weights.push_back(words.Number());
      }
      return spline;
    }

    /**
     *  Reads the words after a curve's kind.
     *
     *  @param  words The text.
     *  @param  kind  The kind.
     *  @return       The curve, or none for an unknown kind.
     */
    auto Read(Reader &words, std::string_view kind) -> std::optional<Curve>
    {
      if (kind == "line")
      {
        const Point start = words.At();
        return Line{.Start = start, .End = words.At()};
      }
      if (kind == "ray")
      {
        const Point start = words.At();
        return Ray{.Start = start, .Way = words.Way()};
      }
      if (kind == "xline")
      {
        const Point through = words.At();
        return Xline{.Through = through, .Way = words.Way()};
      }
      if (kind == "arc")
      {
        const Point  center = words.At();
        const double radius = words.Number();
        const Angle  start  = words.Turn();
        return Arc{.Center = center, .Radius = radius, .Start = start, .Sweep = words.Turn()};
      }
      if (kind == "ellipse")
      {
        const Vector x      = words.Way();
        const Vector y      = words.Way();
        const Point  origin = words.At();
        const Angle  start  = words.Turn();
        return Ellipse{
            .Map = {.X = x, .Y = y, .Origin = origin}, .Start = start, .Sweep = words.Turn()};
      }
      if (kind == "polyline")
      {
        return ReadPolyline(words);
      }
      if (kind == "spline")
      {
        return ReadSpline(words);
      }
      if (kind == "quadratic")
      {
        const Point start   = words.At();
        const Point control = words.At();
        return Bezier::Quadratic{.Start = start, .Control = control, .End = words.At()};
      }
      if (kind == "cubic")
      {
        const Point start = words.At();
        const Point out   = words.At();
        const Point into  = words.At();
        return Bezier::Cubic{.Start = start, .Out = out, .In = into, .End = words.At()};
      }
      return std::nullopt;
    }
  } // namespace

  auto ToText(const Curve &curve) -> std::string
  {
    Writer out;
    std::visit([&out](const auto &shape) -> void { Write(out, shape); }, curve);
    return out.Text();
  }

  auto FromText(std::string_view text) -> std::optional<Curve>
  {
    Reader                 words{text};
    const std::string_view kind  = words.Word();
    std::optional<Curve>   curve = Read(words, kind);
    if (!curve || !words.Done())
    {
      return std::nullopt;
    }
    return curve;
  }
} // namespace Geoxy
