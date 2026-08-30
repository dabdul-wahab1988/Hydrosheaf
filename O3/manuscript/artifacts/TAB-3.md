| Benchmark | Independent / screening condition | Value | Calibrated / prior-informed / emulated condition | Value | Independence lost |
|---|---|---:|---|---:|---|
| Age / residence time | Held-out uncalibrated log10 R2 | 0.482 | Calibrated emulation log10 R2 | 0.962 | Ridge calibration fit on the same held-out folds it is scored against |
| Age / residence time | Held-out uncalibrated within-factor-2 | 0.564 | Calibrated emulation within-factor-2 | 0.769–0.993 | Same non-independence as above |
| Reaction | Chance level, 4-class MRS | 0.250 | Held-out archetype transfer, 4-class MRS accuracy | 0.489 | Trained on 3 archetypes, tested on 1 untouched archetype (a genuine transfer test, not an emulation) |
| Topology | Independent F1 (no MODPATH access) | 0.618 | Prior-informed output-graph edges (override / merge / only) | 302 / 329 / 174 | MODPATH edges entered the graph as prior information; not independent validation |

Notes. The two calibration exercises are not equally rigorous and are not
presented as such: the age benchmark's calibrated value is an emulation fit
and scored on the same data, which is why its improvement (0.482 to 0.962)
is large; the reaction benchmark's Mechanism Resolution Score is tested on a
fourth, untouched archetype it never trained on, a harder and more
independent test, which is why its absolute accuracy (0.489) is far lower
and still only doubles chance. The topology benchmark's prior-informed row
is not a calibration exercise at all -- it is a different mode of use
(MODPATH edges as priors, not as an accuracy target) and is reported on its
own edge-count scale rather than forced onto the same axis.
