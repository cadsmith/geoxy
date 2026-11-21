import { createRequire } from 'node:module'
import { defineConfig } from 'vitepress'

// The pages live outside this folder, in ../docs, where no node_modules is found: vue is named here.
const require = createRequire(import.meta.url)
const vue = (path: string) => require.resolve(path)

const page = (section: string, slug: string, text: string) => ({ text, link: `/${section}/${slug}/` })

const guide = [
  page('guide', 'start', 'Getting started'),
  page('guide', 'install', 'Installing geoxy'),
  page('guide', 'curves', 'Curves'),
  page('guide', 'tolerance', 'Tolerance and accuracy'),
  page('guide', 'errors', 'Errors'),
]

const tutorials = [
  page('tutorials', 'intersecting', 'Where curves meet'),
  page('tutorials', 'snapping', 'Snapping'),
  page('tutorials', 'measuring', 'Measuring'),
  page('tutorials', 'editing', 'Editing curves'),
  page('tutorials', 'offsetting', 'Offsetting'),
  page('tutorials', 'regions', 'Regions'),
  page('tutorials', 'transforming', 'Moving and scaling'),
  page('tutorials', 'indexing', 'Indexing a drawing'),
  page('tutorials', 'displaying', 'Drawing curves'),
]

const concepts = [
  page('concepts', 'exactness', 'Exactness'),
  page('concepts', 'quality', 'How it is tested'),
  page('concepts', 'limits', 'Limits'),
]

export default defineConfig({
  title: 'geoxy',
  description: 'Exact and fast 2D geometry for CAD, in C++23.',
  srcDir: '../docs',
  srcExclude: ['cmake.md', 'commit.md', 'cpp.md'],
  outDir: './.vitepress/dist',
  cacheDir: './.vitepress/cache',
  head: [['link', { rel: 'icon', type: 'image/svg+xml', href: '/favicon.svg' }]],
  vite: {
    publicDir: new URL('../public', import.meta.url).pathname,
    resolve: {
      alias: [
        { find: /^vue\/server-renderer$/, replacement: vue('vue/server-renderer') },
        { find: /^vue$/, replacement: vue('vue') },
      ],
    },
  },
  rewrites: { ':section/:slug.md': ':section/:slug/index.md' },
  themeConfig: {
    logo: '/logo.svg',
    nav: [
      { text: 'Guide', link: '/guide/start/', activeMatch: '/guide/' },
      { text: 'Tutorials', link: '/tutorials/intersecting/', activeMatch: '/tutorials/' },
      { text: 'Concepts', link: '/concepts/exactness/', activeMatch: '/concepts/' },
    ],
    sidebar: [
      { text: 'Guide', items: guide },
      { text: 'Tutorials', items: tutorials },
      { text: 'Concepts', items: concepts },
    ],
    outline: { level: 2 },
    search: { provider: 'local' },
    footer: { message: 'Released under the MIT License.' },
  },
})
