# Eigenbasis

Rebuild a function from eigenmodes, one term at a time. It has three views:

- **Expand.** Watch the partial sums S_N close in on f. The page shows the Gibbs overshoot, Parseval's sum filling up, the size of every coefficient, and how fast they fall.
- **Completeness.** The kernel K_N(x, x₀) = Σ φₙ(x) φₙ(x₀) sharpens into a delta as N grows. Leave one mode out and the error stalls for good.
- **Evolve.** ψ(x, t) = Σ cₙ e^(−iEₙt) φₙ(x): each mode only turns its phase. The view shows the probability density coloured by phase, the quantum carpet over one period, and the chance of finding the starting state again.

There are five bases: Fourier (a particle on a ring), the particle in a box, Legendre (a rigid rotor), the harmonic oscillator, and Haar wavelets, included as a contrast.

## Open it

Double-click `index.html`. Everything it needs is in this folder, so it also works offline, for example on a classroom projector with no Wi-Fi.

## Put it online for free with GitHub Pages

1. Sign in to GitHub and create a new **public** repository, for example `eigenbasis`.
2. In the repository, choose **Add file → Upload files**. Drag in everything inside this folder (`index.html`, `engine.js`, `favicon.svg`, `og-image.png`, `README.md` and the `fonts` folder), then click **Commit changes**.
3. Open **Settings → Pages**. Under **Build and deployment**, set **Source** to **Deploy from a branch**, pick the branch `main` and the folder `/ (root)`, then click **Save**.

After a minute or two the site is live at `https://<your-username>.github.io/eigenbasis/`.

## Link previews on WhatsApp and elsewhere

Preview cards need the full address of the image. Once the site is live, open `index.html` and change

```html
<meta property="og:image" content="og-image.png">
```

to the full URL, for example `https://<your-username>.github.io/eigenbasis/og-image.png`.

## Links to a specific view

The address bar always describes what is on screen, and **Copy link** copies it. Some views worth showing someone:

| Link ending | Shows |
| --- | --- |
| `#/fourier/square/23` | The square wave and its 9% Gibbs overshoot |
| `#/fourier/exp/41` and `#/legendre/exp/16` | The same smooth function: 1/n decay with Fourier (the ends don't match), rounding level after 14 terms with Legendre |
| `#/legendre/runge/60` | Geometric decay at the rate set by the poles at ±i/5 |
| `#/oscillator/shifted/12/p` | A coherent state: its energy probabilities are Poisson with mean 2 |
| `#/haar/box/4` | Exact with 4 terms, and no overshoot ever |
| `#/fourier/square/41/kernel/x0=0.5` | The Dirichlet kernel sitting on the jump |
| `#/legendre/runge/40/kernel/x0=0.985` | The kernel near an end of the interval |
| `#/box/shifted/80/evolve/t=0.25` | A fractional revival in the box: two half-size copies |
| `#/box/gaussian/60/evolve/k0=50` | A kicked packet bouncing between the walls; press play |
| `#/oscillator/shifted/40/evolve` | A coherent state swinging without spreading; press play |
| `#/legendre/f:abs(x)/40` | A typed formula |

The link pattern is `#/basis/function/N`, followed by any options:

- **Bases:** `fourier`, `box`, `legendre`, `oscillator`, `haar`.
- **Functions:** `square`, `box`, `tent`, `ramp`, `gaussian`, `shifted`, `runge`, `exp`, `cusp`, `f:<formula>`. A hand-drawn function is encoded in the link too.
- **Options:** `p` shows coefficient shares, `kernel` or `evolve` picks the view, and `x0=`, `drop=`, `k0=` and `t=` set the marked point, the left-out mode, the push and the time (as a fraction of the period).

## Formulas

Use `x` as the variable. You can write `+ − * / ^`, `|x|`, `π` or `pi`, `e`, implicit products such as `3x` or `2pi x`, and the usual functions (`sin`, `cos`, `exp`, `log`, `sqrt`, `abs`, `sign`, `floor`, `sech`, `min(a, b)`, `max(a, b)` and more). Comparisons give 1 or 0, so `(x>0.3)(x<0.6)` is a box.

## Keyboard

- `←` and `→` add or remove one term; hold `Shift` for ten.
- `Space` plays: terms in Expand and Completeness, time in Evolve.
- `[` and `]` step the time.
- `1` to `5` switch the basis.

## How the numbers are made

Every basis is orthonormal on its interval, so cₙ = ∫ f φₙ dx. The integrals use composite Gauss–Legendre quadrature, with panels split adaptively wherever f jumps or kinks, and are accurate to about 10⁻¹⁴. The page checks this against closed forms: the square wave, the Dirichlet kernel's height N/2π, and a coherent state's Poisson weights and ⟨E⟩ = 2.5. The error ‖f − S_N‖ is integrated directly, not inferred from Parseval, so it stays accurate down to rounding level.

## Credits

STIX Two Text (The STIX Fonts Project) and Atkinson Hyperlegible (Braille Institute of America) are used under the SIL Open Font License 1.1; the `OFL` files are in `fonts/`.
