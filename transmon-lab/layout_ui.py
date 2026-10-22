"""Streamlit UI for process-aware, preliminary GDS and solver handoff."""

from dataclasses import asdict
import json

import streamlit as st

from chip_layout import (
    Process, Readout, compare_extraction, generate_layout, load_process, render_svg,
)


def render_layout_tab(explore_ej: float, explore_ec: float) -> None:
    st.subheader("From a target to a chip layout")
    st.write("Generate a preliminary Nb/Al layout on silicon, inspect it in KLayout, and carry its geometry and ports into HFSS or Q3D.")
    st.info("This export is a starting geometry. The pad/coupler calibrations are examples until replaced with lab data, and the Al junction still needs a process-specific junction cell.")
    options=["Explore circuit"]
    if "design" in st.session_state:
        options.insert(0,"Saved target design")
    source=st.selectbox("Qubit source",options,key="layout_source")
    if source=="Saved target design":
        target_f,target_alpha,result=st.session_state["design"]
        source_ej,source_ec=result.ej_ghz,result.ec_ghz
        st.caption(f"Uses the saved {target_f:g} GHz / {target_alpha:g} MHz target. Readout estimates are for ng = 0.")
    else:
        source_ej,source_ec=explore_ej,explore_ec
        st.caption("Uses EJ and EC from the sidebar, at ng = 0. To lay out a target design, first click Find a design in the Design tab.")

    uploaded=st.file_uploader("Optional: load a process/calibration JSON",type=["json"],key="process_upload")
    if uploaded is not None and st.button("Load process values"):
        try:
            process=load_process(json.loads(uploaded.getvalue()))
            for name,value in asdict(process).items():
                st.session_state["process_"+name]=value
            st.success("Loaded process values. Inspect them below, then generate a new layout.")
        except (ValueError,TypeError) as error:
            st.error(str(error))

    defaults=asdict(Process())

    def number(name,label,minimum,maximum,step):
        key="process_"+name
        if key not in st.session_state:
            st.session_state[key]=defaults[name]
        return st.number_input(label,min_value=minimum,max_value=maximum,step=step,key=key)

    process_values=dict(defaults)
    with st.form("layout_form"):
        a,b,c=st.columns(3)
        frequency=a.number_input("Readout target (GHz)",5.0,10.0,7.0,0.1,key="readout_frequency")
        coupling=b.number_input("Target g/2π (MHz)",5.0,150.0,80.0,5.0,key="readout_g")
        linewidth=c.number_input("Loaded κ/2π (MHz)",0.05,10.0,1.5,0.05,key="readout_kappa")
        a,b=st.columns(2)
        qi=a.number_input("Assumed resonator intrinsic Q",1000.0,10000000.0,200000.0,10000.0)
        drive_cap=b.number_input("Target drive coupling C (fF)",0.01,0.5,0.15,0.01)
        with st.expander("Wafer, stack, and GDS layers",expanded=False):
            a,b,c=st.columns(3)
            with a:
                process_values["chip_width_um"]=number("chip_width_um","Die width (µm; 10000 = 10 mm)",8000.0,20000.0,500.0)
                process_values["chip_height_um"]=number("chip_height_um","Die height (µm)",6000.0,20000.0,500.0)
                process_values["wafer_diameter_mm"]=number("wafer_diameter_mm","Wafer diameter (mm; 50.8 = 2 inch)",25.0,200.0,0.1)
                process_values["wafer_edge_exclusion_mm"]=number("wafer_edge_exclusion_mm","Wafer edge exclusion (mm)",0.1,10.0,0.1)
                process_values["dicing_street_um"]=number("dicing_street_um","Dicing street (µm)",50.0,2000.0,50.0)
            with b:
                process_values["substrate_eps_r"]=number("substrate_eps_r","Silicon relative permittivity",2.0,20.0,0.05)
                process_values["substrate_thickness_um"]=number("substrate_thickness_um","Silicon thickness (µm)",100.0,2000.0,25.0)
                process_values["nb_thickness_nm"]=number("nb_thickness_nm","Nb thickness (nm)",20.0,1000.0,10.0)
                process_values["al_thickness_nm"]=number("al_thickness_nm","Al thickness (nm)",10.0,500.0,5.0)
            with c:
                process_values["nb_layer"]=number("nb_layer","Nb GDS layer",1,80,1)
                process_values["al_layer"]=number("al_layer","Al landing-lead GDS layer",1,80,1)
                process_values["minimum_nb_feature_um"]=number("minimum_nb_feature_um","Minimum Nb feature (µm)",0.1,10.0,0.1)
                process_values["minimum_nb_gap_um"]=number("minimum_nb_gap_um","Minimum Nb gap (µm)",0.1,10.0,0.1)
            st.caption("Thicknesses are exported for EM setup. The analytical wavelength estimate uses a dielectric half-space and omits film/kinetic-inductance corrections. The Al z interval is an editable solver convention; actual junction/contact stacking must follow your process.")
        with st.expander("Geometry calibration — replace example values with lab data",expanded=False):
            key="process_calibration_source"
            if key not in st.session_state:
                st.session_state[key]=defaults["calibration_source"]
            process_values["calibration_source"]=st.text_input("Calibration source / run identifier",key=key)
            a,b=st.columns(2)
            with a:
                process_values["reference_arm_reach_um"]=number("reference_arm_reach_um","Reference cross arm reach (µm)",30.0,1000.0,5.0)
                process_values["reference_arm_width_um"]=number("reference_arm_width_um","Reference cross arm width (µm)",5.0,100.0,1.0)
                process_values["reference_moat_gap_um"]=number("reference_moat_gap_um","Reference ground moat (µm)",2.0,100.0,1.0)
                process_values["reference_pad_capacitance_ff"]=number("reference_pad_capacitance_ff","Reference pad capacitance (fF)",5.0,500.0,1.0)
                process_values["cpw_width_um"]=number("cpw_width_um","CPW center width (µm)",2.0,50.0,0.5)
                process_values["cpw_gap_um"]=number("cpw_gap_um","CPW slot gap (µm)",1.0,30.0,0.5)
            with b:
                process_values["qubit_resonator_gap_um"]=number("qubit_resonator_gap_um","Qubit–resonator gap (µm)",1.0,30.0,0.5)
                process_values["qubit_resonator_cap_ff_per_um"]=number("qubit_resonator_cap_ff_per_um","Qubit–resonator C per overlap (fF/µm)",0.005,1.0,0.005)
                process_values["resonator_feed_gap_um"]=number("resonator_feed_gap_um","Resonator–feedline gap (µm)",1.0,30.0,0.5)
                process_values["resonator_feed_cap_ff_per_um"]=number("resonator_feed_cap_ff_per_um","Resonator–feedline C per overlap (fF/µm)",0.005,1.0,0.005)
                process_values["drive_gap_um"]=number("drive_gap_um","Drive–qubit gap (µm)",1.0,30.0,0.5)
                process_values["drive_cap_ff_per_um"]=number("drive_cap_ff_per_um","Drive C per overlap (fF/µm)",0.001,0.1,0.001)
            a,b=st.columns(2)
            with a:
                process_values["junction_jc_a_cm2"]=number("junction_jc_a_cm2","Al junction critical-current density (A/cm²)",1.0,1000.0,1.0)
            with b:
                process_values["junction_cap_ff_per_um2"]=number("junction_cap_ff_per_um2","Junction specific capacitance (fF/µm²)",1.0,200.0,1.0)
            st.caption("Pad dimensions scale uniformly from the reference cross. Coupling finger lengths use the entered C-per-overlap calibrations at the entered gaps; changing a gap does not automatically recalibrate them. JJ area is Ic/Jc and is exported as an intent marker, not a fabrication junction cell.")
        submitted=st.form_submit_button("Generate chip layout",type="primary")
    if submitted:
        try:
            process=load_process(process_values)
            with st.spinner("Routing, checking net isolation, and building GDS exports…"):
                bundle=generate_layout(source_ej,source_ec,Readout(frequency,coupling,linewidth,qi,drive_cap),process)
            st.session_state["layout_bundle"]=bundle
            st.session_state.pop("layout_error",None)
        except (ValueError,TypeError) as error:
            st.session_state.pop("layout_bundle",None)
            st.session_state["layout_error"]=str(error)
    if "layout_error" in st.session_state:
        st.error(st.session_state["layout_error"])
    if "layout_bundle" not in st.session_state:
        st.caption("Generate a layout to preview its geometry and enable GDS/solver downloads.")
        return

    bundle=st.session_state["layout_bundle"]
    sizing=bundle.geometry["sizing"]
    st.success("Generated layout: no unintended net overlaps, metal inside the die, and bare resonator centerline length checked.")
    st.caption(f"Displayed export: {sizing['target_f01_ghz']:.4f} GHz qubit, {bundle.geometry['readout_targets']['frequency_ghz']:g} GHz bare readout target. Control changes take effect after Generate chip layout. Geometry ID: {bundle.geometry['geometry_sha256'][:16]}…")
    a,b,c,d=st.columns(4)
    a.metric("Bare λ/2 centerline",f"{sizing['bare_half_wave_length_um']/1000:.3f} mm")
    b.metric("CPW impedance estimate",f"{sizing['cpw_impedance_ohm']:.1f} Ω")
    c.metric("Required critical current",f"{sizing['critical_current_na']:.2f} nA")
    d.metric("Complete dies / wafer",str(bundle.manifest["wafer_die_count"]))
    view=st.selectbox("Preview region",["chip","qubit","junction"],key="layout_view")
    st.image(render_svg(bundle.geometry,view),width="stretch")
    st.caption("False colors: gray ground · purple qubit · blue resonator · orange readout feedline · green drive · pink Al landings · red JJ area intent. Geometry is in µm. Ground bridges/bonds and a process junction are not included.")
    a,b,c=st.columns(3)
    a.download_button("Download chip GDS",bundle.chip_gds,"transmon-chip.gds","application/octet-stream")
    b.download_button("Download wafer GDS",bundle.wafer_gds,"transmon-wafer.gds","application/octet-stream")
    c.download_button("Download solver package",bundle.package_zip,"transmon-layout-package.zip","application/zip")
    st.caption("The package contains GDS, DXF construction surfaces, geometry/layer/stack/port JSON, a linear SPICE circuit, and an extraction-results template. HFSS/Q3D still need your conductor, junction, boundary, and excitation assignments; no Ansys solver has been run by this app.")
    with st.expander("Circuit budget and readout estimates"):
        st.dataframe([
            {"Element":"Pad to ground","C (fF)":sizing["pad_capacitance_ff"]},
            {"Element":"Junction","C (fF)":sizing["junction_capacitance_ff"]},
            {"Element":"Qubit–readout","C (fF)":sizing["qubit_resonator_capacitance_ff"]},
            {"Element":"Drive–qubit","C (fF)":sizing["drive_capacitance_ff"]},
            {"Element":"Total qubit diagonal budget","C (fF)":sizing["total_capacitance_ff"]},
        ],hide_index=True,width="stretch")
        st.write(f"Resonator–feedline coupling: {sizing['resonator_feed_capacitance_ff']:.3f} fF. Effective resonator mode capacitance: {sizing['resonator_mode_capacitance_ff']:.1f} fF.")
        if sizing["dispersive_estimate_valid"]:
            st.write(f"Weak-coupling estimates: χ/2π ≈ {sizing['chi_mhz']:.3f} MHz; two-level Purcell-limited T1 ≈ {sizing['purcell_t1_estimate_us']:.1f} µs. No Purcell filter is included.")
        else:
            st.warning("The requested qubit/readout parameters are outside the dispersive criterion used here. χ and Purcell estimates are withheld; coupled-mode analysis is needed.")
        st.write(f"Nominal junction overlap area: {sizing['junction_area_um2']:.4f} µm²; linearized LJ: {sizing['junction_linear_inductance_nh']:.2f} nH.")
        st.caption("The half-wave length is unloaded. Bends, coupling, material inductance, packaging, and extracted effective capacitance change the final frequencies.")
        st.write("Calibration provenance: "+bundle.manifest["calibration_source"])
    st.subheader("Close the loop with extracted values")
    st.write("Fill the package's em_results_template.json after your HFSS/Q3D run. Keep its geometry ID and identify the solver run. The app will compare your extracted values with this layout's targets.")
    extraction=st.file_uploader("Import extraction results JSON",type=["json"],key="em_results_upload")
    if extraction is not None:
        try:
            compared=compare_extraction(bundle,json.loads(extraction.getvalue()))
            st.caption("User-supplied extraction source: "+compared["source"])
            a,b=st.columns(2)
            a.metric("Predicted f01 from extracted C",f"{compared['extracted_f01_ghz']:.4f} GHz",f"{compared['frequency_error_mhz']:+.3f} MHz vs target",delta_color="off")
            b.metric("Predicted α from extracted C",f"{compared['extracted_alpha_mhz']:.2f} MHz",f"{compared['alpha_error_mhz']:+.3f} MHz vs target",delta_color="off")
            st.json(compared)
        except (ValueError,TypeError) as error:
            st.error(str(error))
