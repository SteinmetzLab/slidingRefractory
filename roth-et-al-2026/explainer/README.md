# Sliding RP explainer: movie and interactive explorer

Teaching material for talks, not part of the paper's analyses.

## Outputs (copied to `D:\Dropbox\papers\2026_SlidingRP\explainer\`)

| File | What it is |
|---|---|
| `sliding_rp_explainer_captions.mp4` | Fig 2 built step by step, with a caption for each step. 1920 x 1080, 30 fps, H.264, about 63 s. Drops into PowerPoint. |
| `sliding_rp_explainer_nocaptions.mp4` | The same without captions, to narrate over. |
| `sliding_rp_explorer.html` | Standalone interactive page (one file, no network needed). Sliders for the simulated unit and the metric settings; ACG, confidence curve, single-window Poisson test, confidence matrix; Sliding RP and Hill–Llobet verdicts; 100-train acceptance rates. |

## The movie

`python make_movie.py` (captioned) or `python make_movie.py --no-captions`.
`--stills 40,300,1500` writes PNG stills of chosen frames for checking layout.

It uses the example neuron from the published Fig 2
(`fig2_example_neuron_spike_times.npy`: 10 spikes/s, 2.5 ms refractory period,
plus a 1 spike/s Poisson contaminant, 1 h). Every number shown is computed with
the package's formulas, and the script asserts agreement with
`slidingRP.metrics.computeMatrix` and `slidingRP` before drawing.

Two differences from the printed Fig 2 worth knowing:

* At 1.5 ms the current package counts 110 violations; the figure printed 109
  (75 at 1 ms and 157 at 2.2 ms agree). The figure was drawn with an earlier
  ACG routine.
* The figure's printed probabilities (20.2%, 9.0%, 3.1%) do not follow from its
  own printed counts and means (they give 21.3%, 9.3%, 3.2%). The movie shows
  21.0%, 10.6% and 3.1%, i.e. confidences of 79.0%, 89.4% and 96.9%. So in the
  movie the 1.5 ms example sits just below the 90% line rather than just above.
  Fig 2's annotations should be regenerated from code for the revision.

## The explorer

`srp_core.js` is the numerical core (generator, ACG, Poisson tails, Sliding RP,
C_min, confidence matrix, Hill–Llobet). It is tested against the Python package:

```
python export_reference.py   # writes reference answers to D:/temp/slidingRP_resub/explainer_ref
node test_core.js            # compares; exits non-zero on any failure
python build_explorer.py     # inlines the core into explorer_template.html
```

What the test establishes (last run 2026-09-26, all passing): ACGs identical to
`computeACG`; maximum confidence within 1e-9; C_min within 1e-7 and at the same
bin; the full confidence matrix within 1e-9; first and last accepted windows
identical; Hill–Llobet counts, estimates and verdicts identical; Poisson tails
within 1e-9 up to counts of a million; and acceptance rates of the JavaScript
simulator within sampling error of the Python simulator over 1000 trains.

The page follows the package's conventions: contamination is the fraction of
all spikes from a Poisson source with no refractory period; the Hill–Llobet
window is rounded up to a whole ACG bin, as `RPmetric_Classic.m` does.
