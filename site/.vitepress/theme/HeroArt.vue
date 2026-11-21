<template>
  <!-- A plate with a round corner, its offsets either way, and the place on it nearest a point. -->
  <svg class="hero-art" viewBox="40 60 420 310" role="img"
       aria-label="A shape with its inward and outward offsets and the nearest point to a cursor">
    <defs>
      <pattern id="hero-dots" width="20" height="20" patternUnits="userSpaceOnUse">
        <circle cx="10" cy="10" r="1" class="dot" />
      </pattern>
    </defs>
    <rect x="40" y="60" width="420" height="310" fill="url(#hero-dots)" />
    <path class="band" d="M120 328 H330 A28 28 0 0 0 358 300 V200 A88 88 0 0 0 270 112 H120 A28 28 0 0 0 92 140 V300 A28 28 0 0 0 120 328 Z M148 272 V168 H270 A32 32 0 0 1 302 200 V272 Z" />
    <path class="offset outer" d="M120 328 H330 A28 28 0 0 0 358 300 V200 A88 88 0 0 0 270 112 H120 A28 28 0 0 0 92 140 V300 A28 28 0 0 0 120 328 Z" />
    <path class="offset inner" d="M148 272 H302 V200 A32 32 0 0 0 270 168 H148 Z" />
    <path class="shape" d="M120 300 H330 V200 A60 60 0 0 0 270 140 H120 Z" />
    <g class="vertices">
      <rect x="116" y="296" width="8" height="8" /><rect x="326" y="296" width="8" height="8" />
      <rect x="326" y="196" width="8" height="8" /><rect x="266" y="136" width="8" height="8" />
      <rect x="116" y="136" width="8" height="8" />
    </g>
    <line class="reach" x1="404" y1="86" x2="315.8" y2="161.2" />
    <circle class="foot" cx="315.8" cy="161.2" r="5" />
    <circle class="cursor" cx="404" cy="86" r="4" />
    <text class="label" x="412" y="80" text-anchor="end" dy="-10">Nearest</text>
    <text class="label gold" x="225" y="350">Buffer · Outward</text>
    <text class="label gold" x="225" y="262" text-anchor="middle">Inward</text>
  </svg>
</template>

<style scoped>
.hero-art { width: 100%; max-width: 520px; height: auto; }
.dot { fill: var(--vp-c-divider); }
.band { fill: var(--geoxy-gold); fill-opacity: 0.1; fill-rule: evenodd; }
.shape { fill: none; stroke: var(--vp-c-text-1); stroke-width: 2.5; stroke-linejoin: round; }
.offset { fill: none; stroke: var(--geoxy-accent); stroke-width: 2; stroke-dasharray: 7 5; animation: march 1.6s linear infinite; }
.offset.inner { animation-direction: reverse; }
.vertices rect { fill: var(--vp-c-bg); stroke: var(--vp-c-text-1); stroke-width: 1.6; }
.reach { stroke: var(--vp-c-text-2); stroke-width: 1.4; stroke-dasharray: 3 4; }
.foot { fill: var(--geoxy-accent); stroke: var(--vp-c-bg); stroke-width: 2; }
.cursor { fill: var(--vp-c-text-1); }
.label { font-family: var(--vp-font-family-mono); font-size: 12px; fill: var(--vp-c-text-2); }
.label.gold { fill: var(--geoxy-accent); }
@keyframes march { to { stroke-dashoffset: -24; } }
@media (prefers-reduced-motion: reduce) { .offset { animation: none; } }
</style>
