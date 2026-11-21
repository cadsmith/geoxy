# C++ Guidelines

## Naming

- Namespaces, types, functions, members and constants are PascalCase.
- Parameters and locals are camelCase words; a single letter only for math.
- Private members start with `m_`.
- Template type parameters end in `T`.
- Binary operators name their operands `lhs` and `rhs`.

```cpp
namespace Foo
{
  template <typename BarT>
  class Baz
  {
  public:
    auto Qux(int quux) -> int;

  private:
    int m_Corge = 0;
  };

  auto operator+(Baz<int> lhs, Baz<int> rhs) -> Baz<int>;
} // namespace Foo
```

## Files

- A file holds one responsibility and is named by one lowercase word.
- A directory replaces a compound name; namespaces follow directories.

```txt
foo/bar.hpp      Foo::Bar
foo/baz.hpp      Foo::Baz
foo/bar/qux.hpp  Foo::Bar::Qux
```

## Functions and types

- A function takes at most three parameters.
- A class has one responsibility.
- A comparison is a hidden friend.

```cpp
struct Foo
{
  int Bar = 0;
  int Baz = 0;

  friend auto operator==(const Foo &lhs, const Foo &rhs) -> bool;
};
```

## Values

- A value is initialized with `=`; braces alone call a constructor.
- Floating point values compare through a tolerance, never with `==`.

```cpp
const Foo    foo = {.Bar = 1, .Baz = 2};
const double bar = 0.0;
const Qux    qux{1};
```

## Layout

- A line ends by the 100th column.
- A namespace indents its body.
- Consecutive declarations align on their names and their `=`.
- Whole numbers in a run right-align.
- A list of struct rows takes one row a line, its fields in columns.
- A broken condition takes one operand a line, the operator leading, aligned on
  its comparisons.
- Consecutive assertions align on their `==`.

```cpp
namespace Foo
{
  auto Bar(const Baz &baz, Qux qux) -> bool
  {
    const int    corge  =     -1;
    const int    grault = 20'000;
    const double garply = 0.5;

    return baz.Low.Fred <= qux.Fred
        && qux.Fred     <= baz.High.Fred
        && baz.Low.Thud <= qux.Thud
        && qux.Thud     <= baz.High.Thud;
  }

  constexpr std::array Waldos{
      Waldo{.Fred =    0, .Thud =     1},
      Waldo{.Fred = -5e6, .Thud = 1e-20},
      Waldo{.Fred =  1e8, .Thud =     2},
  };
} // namespace Foo
```

```cpp
CHECK(Foo(1)    == Bar);
CHECK(Foo(22)   == Baz);
CHECK(Foo(-333) == Qux);
```

## Comments

- A file opens with a framed banner: `@file`, `@brief`, `@copyright`.
- A declaration takes a doc block: a summary, details after an empty line, then
  its tags as one aligned table.
- Every parameter has an `@param`, every value a `@return`.
- A member or an enumerator takes a trailing `/**< */` note.

```cpp
/***************************************************************************************************
 *  @file      foo/bar.hpp
 *  @brief     Bars.
 *  @copyright Copyright (c) 2026 Foo. Released under the MIT License.
 **************************************************************************************************/

namespace Foo
{
  /**
   *  A bar.
   *
   *  A baz and a qux, held together.
   */
  struct Bar
  {
    int Baz     = 0; /**< The baz.             */
    int QuxQuux = 0; /**< The qux of the quux. */
  };

  /**
   *  The corge of a bar.
   *
   *  @pre           The bar is not empty.
   *  @param  bar    The bar.
   *  @param  grault The grault.
   *  @return        The corge.
   */
  auto Corge(Bar bar, int grault) -> int;

  enum class Fred
  {
    Thud  = -1, /**< The thud.  */
    Plugh =  0, /**< The plugh. */
    Xyzzy = 22, /**< The xyzzy. */
  };
} // namespace Foo
```

## Tooling

Nothing above is kept by hand. clang-format lays out the code, then
`tools/format/style.py` frames, wraps and aligns what clang-format cannot.
Doxygen fails on an undocumented declaration, and `tools/format/editor.py` on
anything clangd would mark. VS Code runs the formatter on save.

```sh
python3 -I tools/format/style.py
python3 -I tools/format/style.py --check
python3 -I tools/format/editor.py
cmake --build --preset debug --target geoxy_docs
```
