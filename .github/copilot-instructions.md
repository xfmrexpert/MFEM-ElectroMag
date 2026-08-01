# Copilot Instructions

## Project Guidelines
- Prefer explicit, intention-revealing computational implementations over reusing unrelated intermediate calculations as shortcuts, even when reuse would reduce coding effort.
- Prefer abstractions only when they materially simplify understanding; avoid both bespoke over-abstraction and low-level spaghetti, and value alignment with conventional MFEM solver structure where sensible (don't be a slave to MFEM per se).
- Ensure that electrostatic, magnetostatic, and magnetoquasistatic solvers retain similar high-level structures to facilitate the extraction of genuinely common abstractions across all three.

## Code Style
- Use ASCII-only characters in code and comments; do not introduce non-ASCII characters.