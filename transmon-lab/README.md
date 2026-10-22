# Transmon Lab

A superconducting-circuit EDA app for Gal: explore a transmon, fit a circuit to target frequency/anharmonicity, and generate a preliminary chip/wafer layout with solver handoff.

## Run

The cloud environment keeps dependencies outside the existing repository:

```bash
cd /workspace/eda-lab
/workspace/.cloud-env/eda-lab/venv/bin/python -m streamlit run app.py
```

The service listens on port 8501. Cloud onboarding does not provide a web preview; use a workspace with browser forwarding to interact with it, or run it on your own machine:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m streamlit run app.py
```

Requires Python 3.12. The dependency lock was tested on Linux x86_64. No API keys or external services are needed.

On Windows, from the app folder in PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.lock
.\.venv\Scripts\python.exe -m streamlit run app.py --server.headless false --server.address 127.0.0.1
```

Alternatively, double-click **Start Transmon Lab.bat**. It creates or reuses the local virtual environment, installs the pinned dependencies, and opens the app in your browser. Python 3.12 must be installed. The launcher has been reviewed but not executed on Windows in this cloud task.

For an existing installation, stop the server with Ctrl+C, update the source files, rerun the pip command, then restart Streamlit. Keep the existing `.venv`; no PowerShell activation or execution-policy change is needed.

## First experiment

1. In **Explore**, increase the total capacitance. Observe how EC, f01, and anharmonicity change.
2. Compare the calculated spectrum with the transmon approximation.
3. In **Design**, try a 5 GHz frequency and −250 MHz anharmonicity. Download the candidate and inspect its numerical checks.
4. Read **Learn the agent workflow**. The optimizer is a numerical tool; a future AI agent can call it and explain its choices.

## Chip layout and GDS

1. Find a target circuit in **Design**.
2. Open **Chip layout & GDS**. The saved design is selected as the qubit source.
3. Set readout frequency, coupling g/2π, linewidth κ/2π, and drive coupling.
4. Inspect the process assumptions. Defaults are a 2-inch silicon wafer, 10 × 10 mm dies, Nb large structures on GDS layer 1, and Al landing leads on layer 2. Substrate/film thicknesses and all calibrations are editable examples.
5. Click **Generate chip layout**. Inspect the chip, qubit and junction views.
6. Download the chip GDS, wafer GDS, or solver package using the app's own download buttons. Open GDS in KLayout.

The first unit contains a grounded cross transmon, an open-open half-wave readout resonator, a through feedline with two bond pads, and a separate capacitive drive line with its bond pad. The selected wafer margins/dicing streets yield 12 complete dies at the defaults; wafer flat/notch and test structures are not modeled.

The generator checks net isolation, die bounds, and bare resonator centerline length. It does **not** run a fabrication DRC or an EM solve. The Al landing leads stop at a nominal junction site. GDS layer 91 is an **area-intent marker**, not a Dolan/Manhattan junction or an executable fabrication recipe. Add your process junction/contact cell before fabrication.

Layers 90–94 are annotation/intent layers, not mask layers. Materials are separate GDS layers; colors in the preview identify circuit nets and do not represent mask layers.

### HFSS and Q3D

The solver package includes:

- Positive-metal chip/wafer GDS and a false-color SVG preview.
- `construction.dxf`: micrometre construction surfaces. Build ground as `GROUND_BLANK` minus the union of `GROUND_SUBTRACT`; other NB/AL layers are conductor surfaces.
- `geometry.json`: geometry hash, polygons, clearance construction geometry, units, stack, material identities, nets, port coordinates and analytical sizing assumptions.
- `ports.csv`, `manifest.json`, `process.json`, and wafer placements.
- `circuit.cir`: a **linearized** junction/shunt circuit with a distributed lossless resonator; it does not predict quantum anharmonicity.
- `em_results_template.json`: enter extracted effective qubit capacitance and optional readout/g/linewidth results with a solver/run identifier. Import it in the app; the geometry hash prevents mixing results from different layouts.

HFSS/Q3D still require conductor, junction, excitation, reference-ground, boundary, and material assignments. There is no Ansys execution in this app. Nb/Al material names and z intervals are handoff inputs, not validated superconducting-loss models. Use your process contact stack and a suitable PEC/surface-impedance/superconductor model.

The analytic CPW model assumes a zero-thickness dielectric half-space; wavelength excludes bend/end/coupler loading and kinetic inductance. Cross-pad capacitance scales from an entered reference. Coupling capacitances use entered per-overlap calibrations at the entered gaps. Default calibration values are examples and do not establish a physical target match. Extract the coupled capacitance matrix/effective EC and modes to close the design loop.

## Model

`transmon.py` diagonalizes `H = 4 EC (n − ng)² − EJ cos(phi)` in a charge basis. EJ/h and EC/h are in GHz; capacitance is in fF; offset charge is in Cooper-pair units. The displayed anharmonicity is `(E2 − 2 E1 + E0)/h`, in MHz.

The explorer uses 51 charge states and compares transition frequencies with 71 states. The design solver uses 71 states, verifies its target residuals, and checks convergence against 91 states. This is an isolated single-junction model: it omits loss, coupling, extra modes, fabrication uncertainty, and coherence predictions.

## Validate

```bash
cd /workspace/eda-lab
/workspace/.cloud-env/eda-lab/venv/bin/python -m pytest -q
```

With the server running, check its actual HTTP and WebSocket session:

```bash
/workspace/.cloud-env/eda-lab/venv/bin/python tests/check_live.py
```

The original scattering scripts remain unchanged. Their latest simulation runs with the base scientific Python runtime; `simulation_v3.py` has an existing scalar/array error at line 155.
