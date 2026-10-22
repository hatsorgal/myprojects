# Changes

## 0.2 — Chip layout and solver handoff

- Generate a preliminary cross-transmon/readout unit on a 10 × 10 mm silicon die: Nb pads/CPW, Al junction landing leads, a half-wave readout resonator, a through feedline, a drive line, and three bond pads.
- Preview chip, qubit, and junction regions; export chip/wafer GDS and a solver package containing DXF construction surfaces, geometry/stack/net/port JSON, and a linear small-signal circuit.
- Pack complete dies on a 2-inch wafer with configurable margins and dicing streets.
- Import process/calibration profiles and back-annotate extracted effective capacitance/readout values using a matching geometry ID.
- Check net isolation, die bounds, bare centerline length, file integrity, and wafer fit. Independently read example exports with KLayout.
- Keep geometry calibration provenance explicit. Defaults are examples, the junction cell is process-specific and absent, and HFSS/Q3D/DRC have not been run.
- Add a Windows double-click launcher. The app still uses numerical tools; no AI provider is connected.

## 0.1 — Transmon explorer

- Charge-basis spectrum, approximate formulas, target-design optimization, basis convergence checks, and an introduction to agent tool use.
