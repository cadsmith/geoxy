# geoxy site

The documentation website: VitePress, built to static files from the pages
in `../docs`, in the brand geoxy shares with dxfio and dwgio: a gold tile,
ink, Inter and JetBrains Mono.

```sh
npm install
npm run dev        # http://localhost:5173
npm run build      # .vitepress/dist, every page a static file
npm run preview    # serves the build
```

## Writing the docs

The guide, the tutorials and the concepts are Markdown, one file a page, in
the repository's `docs/guide/`, `docs/tutorials/` and `docs/concepts/`;
`.vitepress/config.ts` lists them in the order the sidebar shows them.
`docs/index.md` is the home page. In a page, `#` is its title, the text under
it its lead, and each `##` a section of the outline, its id given as
`## Title {#id}`. Fences between `::: code-group` and `:::` become tabs, each
named in brackets: ` ```cmake [CMakeLists.txt] `. A page links another by its
address, such as `[Curves](/guide/curves/)`.

A fence opened as ` ```cpp [main.cpp] ` is a whole program, and the fence
after it, ` ```text [Output] `, what it prints:
`python3 -I tools/docs/examples.py` builds and runs each and fails where the
two differ.

## Brand

`public/logo.svg` is the mark: the gold tile with two axes, their origin and a
circle. `.vitepress/theme/custom.css` sets the colors and the fonts, which
are in `public/fonts/` under the SIL Open Font License, and styles the home
page; `.vitepress/theme/HeroArt.vue` is the drawing beside its title.

## Publishing

Each page has one address, such as `/tutorials/offsetting/`, and
`npm run build` writes it to its own file,
`.vitepress/dist/tutorials/offsetting/index.html`. Any static host serves
that folder.
