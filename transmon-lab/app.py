"""Interactive single-junction transmon explorer for learning EDA and agents."""

from dataclasses import asdict
import json

import altair as alt
import numpy as np
import streamlit as st

from transmon import (
    approximation, charge_dispersion_khz, charging_energy_ghz,
    convergence_error_khz, design_transmon, observables, spectrum,
)
from layout_ui import render_layout_tab


st.set_page_config(page_title="Transmon Lab", page_icon="⚛", layout="wide")
st.markdown("""
<style>
.block-container {max-width: 1250px; padding-top: 2rem;}
h1 {letter-spacing: -0.04em;}
[data-testid="stMetric"] {border: 1px solid #e0e5ee; border-radius: 12px; padding: 16px; background: #fff;}
.eyebrow {color:#465fc7; font-size:0.8rem; letter-spacing:0.15em; font-weight:700;}
.intro {font-size:1.12rem; color:#566176; max-width:740px; margin-bottom:1.6rem;}
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="eyebrow">TRANSMON LAB · GAL’S EDA WORKBENCH</div>', unsafe_allow_html=True)
st.title("From circuit to qubit.")
st.markdown('<div class="intro">Explore a single-junction transmon, understand its spectrum, and work backward from the qubit you want to build.</div>', unsafe_allow_html=True)

with st.sidebar:
    st.header("Your circuit")
    st.caption("A Josephson junction shunted by the total island capacitance CΣ.")
    capacitance = st.slider("Total capacitance CΣ (fF)", 20.0, 200.0, 80.0, 1.0, key="capacitance_ff")
    ej = st.slider("Josephson energy EJ/h (GHz)", 1.0, 60.0, 20.0, 0.5, key="ej_ghz")
    ng = st.slider("Offset charge ng (in units of 2e)", 0.0, 0.5, 0.0, 0.01, key="offset_charge")
    ec = charging_energy_ghz(capacitance)
    st.divider()
    st.metric("Charging energy EC/h", f"{ec:.4f} GHz")
    st.caption("CΣ includes junction, shunt, and other island capacitances. EJ and EC here are energies divided by h.")


@st.cache_data
def sweep(ec_value: float, ng_value: float) -> list[dict]:
    rows = []
    for value in np.linspace(1, 60, 120):
        calculated = observables(float(value), ec_value, ng_value)
        estimated = approximation(float(value), ec_value)
        rows.extend([
            {"EJ/h (GHz)": float(value), "f01 (GHz)": calculated.f01_ghz, "Method": "Charge-basis calculation"},
            {"EJ/h (GHz)": float(value), "f01 (GHz)": estimated.f01_ghz, "Method": "Transmon approximation"},
        ])
    return rows


current = observables(ej, ec, ng)
estimated = approximation(ej, ec)
energies = spectrum(ej, ec, ng)
basis_error = convergence_error_khz(ej, ec, ng)

explore, design, layout, learn = st.tabs(["Explore the circuit", "Design to a target", "Chip layout & GDS", "Learn the agent workflow"])

with explore:
    st.write("")
    columns = st.columns(4)
    columns[0].metric("Qubit frequency · f01", f"{current.f01_ghz:.3f} GHz")
    columns[1].metric("Anharmonicity · α", f"{current.alpha_mhz:.1f} MHz")
    columns[2].metric("Energy ratio · EJ/EC", f"{ej/ec:.1f}")
    dispersion = charge_dispersion_khz(ej, ec)
    columns[3].metric("Charge dispersion · f01", "<0.01 kHz" if dispersion < 0.01 else f"{dispersion:,.2f} kHz")
    st.caption("α = f12 − f01. Charge dispersion compares ng = 0 with ng = 0.5; it does not predict a coherence time.")
    if ej / ec < 20:
        st.warning("This circuit has EJ/EC below 20. Charge sensitivity can be substantial and the transmon approximation is less accurate; use the calculated spectrum.")
    if basis_error > 1:
        st.error("The charge-basis convergence check exceeds 1 kHz. Treat this spectrum as unresolved.")

    left, right = st.columns([1, 1.45], gap="large")
    with left:
        st.subheader("A ladder with unequal steps")
        level_rows = [{"start": 0, "end": 1, "Energy / h (GHz)": float(value),
                       "State": f"|{index}⟩"} for index, value in enumerate(energies)]
        level_chart = alt.Chart(alt.Data(values=level_rows)).mark_rule(strokeWidth=4).encode(
            x=alt.X("start:Q", axis=None), x2="end:Q",
            y=alt.Y("Energy / h (GHz):Q", scale=alt.Scale(zero=True)),
            color=alt.Color("State:N", scale=alt.Scale(range=["#4560c9", "#7c5cca", "#b57ac8", "#d5a3bc"]), legend=None),
            tooltip=["State:N", alt.Tooltip("Energy / h (GHz):Q", format=".4f")],
        )
        labels = alt.Chart(alt.Data(values=level_rows)).mark_text(align="right", dx=-7, dy=-10).encode(
            x=alt.X("end:Q", axis=None), y="Energy / h (GHz):Q", text="State:N")
        st.altair_chart((level_chart + labels).properties(height=300), width="stretch")
        st.caption(f"f01 = {current.f01_ghz:.4f} GHz · f12 = {current.f12_ghz:.4f} GHz. Energies are relative to |0⟩.")
    with right:
        st.subheader("What changes when EJ changes?")
        curve = alt.Chart(alt.Data(values=sweep(ec, ng))).mark_line(strokeWidth=2.5).encode(
            x=alt.X("EJ/h (GHz):Q"), y=alt.Y("f01 (GHz):Q", scale=alt.Scale(zero=False)),
            color=alt.Color("Method:N", scale=alt.Scale(range=["#4560c9", "#d29153"]), legend=alt.Legend(orient="bottom", title=None)),
            strokeDash=alt.StrokeDash("Method:N", legend=None),
            tooltip=[alt.Tooltip("EJ/h (GHz):Q", format=".2f"), alt.Tooltip("f01 (GHz):Q", format=".4f"), "Method:N"],
        )
        marker = alt.Chart(alt.Data(values=[{"EJ/h (GHz)": ej}])).mark_rule(color="#8390a5", strokeDash=[3, 3]).encode(x="EJ/h (GHz):Q")
        st.altair_chart((curve + marker).properties(height=300).interactive(), width="stretch")
        st.caption("Capacitance and offset charge stay fixed. Hover for values; scroll to zoom.")

    with st.expander("The model, approximations, and numerical check"):
        st.latex(r"H = 4E_C(\hat n-n_g)^2 - E_J\cos\hat\varphi, \qquad E_C=\frac{e^2}{2C_\Sigma}")
        st.latex(r"f_{01}\approx\frac{\sqrt{8E_J E_C}-E_C}{h},\qquad\alpha\approx-\frac{E_C}{h}")
        st.write(f"The approximation gives f01 = {estimated.f01_ghz:.4f} GHz and α = {estimated.alpha_mhz:.1f} MHz. Its frequency differs from the calculation by {abs(estimated.f01_ghz-current.f01_ghz)*1000:.2f} MHz.")
        st.write(f"The calculation uses 51 charge states. Increasing to 71 changes the first three transition frequencies by at most {basis_error:.5f} kHz.")
        st.caption("This isolated, single-mode model omits resonator coupling, loss, quasiparticles, fabrication variation, and additional circuit modes. It is an educational design estimate, not a device validation.")

with design:
    st.subheader("Start with the qubit you want")
    st.write("Choose f01 and α. The solver finds EJ and total capacitance using the same Hamiltonian, at ng = 0.")
    with st.form("design_target"):
        c1, c2 = st.columns(2)
        target_f = c1.number_input("Target f01 (GHz)", min_value=2.0, max_value=10.0, value=5.0, step=0.1, key="target_f01")
        target_alpha = c2.number_input("Target α (MHz)", min_value=-500.0, max_value=-50.0, value=-250.0, step=10.0, key="target_alpha")
        submitted = st.form_submit_button("Find a design", type="primary")
    if submitted:
        try:
            with st.spinner("Solving and checking the spectrum…"):
                result = design_transmon(target_f, target_alpha)
            st.session_state["design"] = (target_f, target_alpha, result)
            st.session_state["layout_source"] = "Saved target design"
            st.session_state.pop("design_error", None)
        except ValueError as error:
            st.session_state.pop("design", None)
            st.session_state["design_error"] = str(error)
    if "design_error" in st.session_state:
        st.error(st.session_state["design_error"])
    if "design" in st.session_state:
        saved_f, saved_alpha, result = st.session_state["design"]
        st.success(f"Verified in the model: f01 = {result.predicted.f01_ghz:.4f} GHz and α = {result.predicted.alpha_mhz:.2f} MHz.")
        st.caption(f"Saved result for target {saved_f:g} GHz / {saved_alpha:g} MHz at ng = 0. Sidebar controls affect the Explore tab; this result stays tied to its target.")
        a, b, c = st.columns(3)
        a.metric("Suggested EJ/h", f"{result.ej_ghz:.4f} GHz")
        b.metric("Suggested CΣ", f"{result.capacitance_ff:.2f} fF")
        c.metric("EJ/EC", f"{result.ej_ghz/result.ec_ghz:.1f}")
        if result.ej_ghz / result.ec_ghz < 20:
            st.warning("The fitted circuit has EJ/EC below 20. Consider its charge sensitivity before treating it as a transmon design.")
        st.write(f"Verification: |frequency error| = {abs(result.frequency_error_khz):.4f} kHz; |α error| = {abs(result.alpha_error_khz):.4f} kHz; basis change = {result.basis_error_khz:.5f} kHz.")
        st.caption(f"{result.simulator_calls} simulator calls. These errors measure the fit within this model, not the accuracy of a fabricated device.")
        payload = {"model": "single-junction transmon", "offset_charge_2e": 0,
                   "target_f01_ghz": saved_f, "target_alpha_mhz": saved_alpha, **asdict(result)}
        st.download_button("Download design JSON", json.dumps(payload, indent=2), "transmon-design.json", "application/json")
    else:
        st.info("Try the default 5 GHz, −250 MHz target to produce your first candidate circuit.")

with layout:
    render_layout_tab(ej,ec)

with learn:
    st.subheader("An agent needs a reliable tool")
    st.write("You are building that tool first: a simulator whose inputs, units, outputs, and checks are explicit.")
    st.markdown("""
1. **Set a goal.** For example: “Find a 5 GHz transmon with −250 MHz anharmonicity.”
2. **Call a tool.** An AI agent could call `design_transmon(5.0, -250.0)` instead of inventing circuit values.
3. **Check the result.** Recalculate the spectrum, inspect target errors, and increase the charge-basis size.
4. **Make a decision.** Accept the modeled candidate, change constraints, or ask for more evidence.
""")
    st.code('design_transmon(target_f01_ghz=5.0, target_alpha_mhz=-250.0)', language="python")
    st.info("The current app uses a numerical optimizer. No AI model or autonomous agent is connected, and no API key is needed.")
    st.write("**Your first experiment:** increase CΣ in Explore. Predict what happens to EC, f01, and α, then check the plots. Next, use Design to find a circuit at your own target.")
    with st.expander("Where could reinforcement learning fit?"):
        st.write("Reinforcement learning learns a policy through repeated actions and rewards. A later project might tune a simulated calibration experiment, with actions such as changing a pulse parameter and a reward based on fidelity. For this two-parameter design problem, a numerical solver is a practical starting point.")
    st.write("**The layout loop:** target → circuit model → GDS candidate → HFSS/Q3D extraction → compare against the target → revise. An agent can coordinate these tools, while process calibration and solver results supply the evidence.")
    st.caption("Next agent step: let an AI agent call the design and layout tools, propose a change, and explain the numerical and extracted checks before accepting it.")
