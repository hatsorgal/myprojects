"""Solver-neutral, preliminary transmon/readout layout generation.

All drawing coordinates are micrometres. Calibration values are explicit inputs;
the generator does not claim that analytical circuit targets are EM-verified.
"""

from dataclasses import asdict, dataclass, fields
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED

import ezdxf
import gdstk
import numpy as np
from scipy.constants import c, elementary_charge, h, hbar, physical_constants
from scipy.special import ellipk

from transmon import capacitance_ff, charging_energy_ghz, observables


PHI0 = physical_constants["mag. flux quantum"][0]


@dataclass(frozen=True)
class Process:
    substrate_eps_r: float = 11.45
    substrate_thickness_um: float = 500.0
    nb_thickness_nm: float = 200.0
    al_thickness_nm: float = 50.0
    nb_layer: int = 1
    al_layer: int = 2
    chip_width_um: float = 10000.0
    chip_height_um: float = 10000.0
    wafer_diameter_mm: float = 50.8
    wafer_edge_exclusion_mm: float = 2.0
    dicing_street_um: float = 500.0
    cpw_width_um: float = 10.0
    cpw_gap_um: float = 6.0
    reference_arm_reach_um: float = 140.0
    reference_arm_width_um: float = 24.0
    reference_moat_gap_um: float = 12.0
    reference_pad_capacitance_ff: float = 80.0
    qubit_resonator_gap_um: float = 10.0
    qubit_resonator_cap_ff_per_um: float = 0.08
    resonator_feed_gap_um: float = 8.0
    resonator_feed_cap_ff_per_um: float = 0.05
    drive_gap_um: float = 10.0
    drive_cap_ff_per_um: float = 0.005
    junction_jc_a_cm2: float = 50.0
    junction_cap_ff_per_um2: float = 50.0
    minimum_nb_feature_um: float = 1.0
    minimum_nb_gap_um: float = 0.5
    calibration_source: str = "EXAMPLE VALUES — replace with lab extraction or measurements"


@dataclass(frozen=True)
class Readout:
    frequency_ghz: float = 7.0
    coupling_g_mhz: float = 80.0
    loaded_linewidth_mhz: float = 1.5
    intrinsic_q: float = 200000.0
    drive_capacitance_ff: float = 0.15


def _positive_fields(instance) -> None:
    for field in fields(instance):
        value = getattr(instance, field.name)
        if field.name == "calibration_source":
            if not isinstance(value, str) or not value.strip():
                raise ValueError("Calibration source must be a nonempty description.")
        elif isinstance(value,bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{field.name} must be finite and positive.")


def load_process(data: dict) -> Process:
    if not isinstance(data, dict):
        raise ValueError("The process file must contain a JSON object.")
    known = {field.name for field in fields(Process)}
    unknown = set(data) - known
    if unknown:
        raise ValueError("Unknown process fields: " + ", ".join(sorted(unknown)))
    process = Process(**data)
    validate_process(process)
    return process


def validate_process(process: Process) -> None:
    _positive_fields(process)
    if not 2 <= process.substrate_eps_r <= 20:
        raise ValueError("Relative substrate permittivity must be between 2 and 20.")
    if process.chip_width_um < 8000 or process.chip_height_um < 6000:
        raise ValueError("This routing template needs at least an 8 mm × 6 mm die.")
    for layer in [process.nb_layer, process.al_layer]:
        if not isinstance(layer, int) or not 1 <= layer <= 80:
            raise ValueError("Nb and Al layer numbers must be integers from 1 to 80.")
    if process.nb_layer == process.al_layer:
        raise ValueError("Nb and Al must use distinct layers.")
    if process.reference_arm_width_um >= process.reference_arm_reach_um:
        raise ValueError("The reference cross must have arms longer than their width.")
    if min(process.cpw_width_um, process.reference_arm_width_um) < process.minimum_nb_feature_um:
        raise ValueError("A Nb conductor is narrower than the specified minimum feature.")
    if min(process.cpw_gap_um, process.qubit_resonator_gap_um,
           process.resonator_feed_gap_um, process.drive_gap_um) < process.minimum_nb_gap_um:
        raise ValueError("A layout gap is below the specified Nb minimum gap.")
    ranges={
        "substrate_thickness_um":(100,2000),"nb_thickness_nm":(20,1000),"al_thickness_nm":(10,500),
        "chip_width_um":(8000,20000),"chip_height_um":(6000,20000),
        "wafer_diameter_mm":(25,200),"wafer_edge_exclusion_mm":(0.1,10),"dicing_street_um":(50,2000),
        "cpw_width_um":(2,50),"cpw_gap_um":(1,30),
        "reference_arm_reach_um":(30,1000),"reference_arm_width_um":(5,100),
        "reference_moat_gap_um":(2,100),"reference_pad_capacitance_ff":(5,500),
        "qubit_resonator_gap_um":(1,30),"qubit_resonator_cap_ff_per_um":(0.005,1),
        "resonator_feed_gap_um":(1,30),"resonator_feed_cap_ff_per_um":(0.005,1),
        "drive_gap_um":(1,30),"drive_cap_ff_per_um":(0.001,0.1),
        "junction_jc_a_cm2":(1,1000),"junction_cap_ff_per_um2":(1,200),
        "minimum_nb_feature_um":(0.1,10),"minimum_nb_gap_um":(0.1,10),
    }
    for name,(lower,upper) in ranges.items():
        if not lower<=getattr(process,name)<=upper:
            raise ValueError(f"{name} must be between {lower} and {upper} for this template.")


def cpw_parameters(process: Process) -> tuple[float, float]:
    """Zero-thickness CPW on an infinite dielectric half-space: (eps_eff, Z0)."""
    eps_eff = (process.substrate_eps_r + 1) / 2
    k = process.cpw_width_um / (process.cpw_width_um + 2 * process.cpw_gap_um)
    impedance = 30 * math.pi / math.sqrt(eps_eff) * ellipk(1 - k*k) / ellipk(k*k)
    return eps_eff, float(impedance)


def size_circuit(ej_ghz: float, ec_ghz: float, readout: Readout, process: Process) -> dict:
    validate_process(process)
    _positive_fields(readout)
    if not 5 <= readout.frequency_ghz <= 10:
        raise ValueError("The layout template supports readout targets from 5 to 10 GHz.")
    qubit = observables(ej_ghz, ec_ghz, charge_cutoff=35)
    total_c = capacitance_ff(ec_ghz)
    eps_eff, z0 = cpw_parameters(process)
    fr = readout.frequency_ghz * 1e9
    omega = 2 * math.pi * fr
    # Open-open half-wave resonator, both coupling sites at voltage antinodes.
    mode_c = math.pi / (2 * omega * z0)
    voltage_zpf = math.sqrt(hbar * omega / (2 * mode_c))
    n01 = (ej_ghz / (32 * ec_ghz)) ** 0.25
    c_qr = readout.coupling_g_mhz * 1e6 * h * (total_c * 1e-15) / (2 * elementary_charge * voltage_zpf * n01) / 1e-15
    inverse_q_external = readout.loaded_linewidth_mhz * 1e6 / fr - 1 / readout.intrinsic_q
    if inverse_q_external <= 0:
        raise ValueError("The requested loaded linewidth is below the intrinsic linewidth f/Qi.")
    q_external = 1 / inverse_q_external
    # A matched through feedline presents Z0/2 at its capacitive connection.
    c_rf = math.sqrt(mode_c / (omega * q_external * (z0 / 2))) / 1e-15
    ic_a = 2 * math.pi * (ej_ghz * 1e9 * h) / PHI0
    junction_area = ic_a / (process.junction_jc_a_cm2 * 1e-8)
    c_junction = process.junction_cap_ff_per_um2 * junction_area
    c_pad = total_c - c_qr - readout.drive_capacitance_ff - c_junction
    if c_pad <= 0:
        raise ValueError("Coupling and junction capacitance exhaust the total qubit capacitance budget.")
    scale = c_pad / process.reference_pad_capacitance_ff
    arm = process.reference_arm_reach_um * scale
    width = process.reference_arm_width_um * scale
    moat = process.reference_moat_gap_um * scale
    if width < process.minimum_nb_feature_um or moat < process.minimum_nb_gap_um:
        raise ValueError("The scaled transmon violates the entered feature/gap limits.")
    qr_overlap = c_qr / process.qubit_resonator_cap_ff_per_um
    rf_overlap = c_rf / process.resonator_feed_cap_ff_per_um
    drive_overlap = readout.drive_capacitance_ff / process.drive_cap_ff_per_um
    usable_arm = arm - width / 2 - process.minimum_nb_gap_um
    if max(qr_overlap, drive_overlap) > usable_arm:
        raise ValueError("A coupling finger does not fit its transmon arm. Reduce the coupling target or supply a calibrated geometry with greater capacitance per length.")
    if rf_overlap > 1500:
        raise ValueError("The feedline coupling section exceeds this template's 1500 µm limit.")
    detuning = qubit.f01_ghz - readout.frequency_ghz
    alpha_ghz = qubit.alpha_mhz / 1000
    g_ghz = readout.coupling_g_mhz / 1000
    dispersive = abs(detuning) > 10*g_ghz and abs(detuning + alpha_ghz) > 10*g_ghz
    chi_mhz = g_ghz**2 * alpha_ghz / (detuning * (detuning + alpha_ghz)) * 1000 if dispersive else None
    purcell_us = (detuning/g_ghz)**2 / (2*math.pi*readout.loaded_linewidth_mhz) if dispersive else None
    return {
        "ej_ghz": ej_ghz, "ec_ghz": ec_ghz, "target_f01_ghz": qubit.f01_ghz,
        "target_alpha_mhz": qubit.alpha_mhz, "total_capacitance_ff": total_c,
        "pad_capacitance_ff": c_pad, "junction_capacitance_ff": c_junction,
        "qubit_resonator_capacitance_ff": c_qr, "resonator_feed_capacitance_ff": c_rf,
        "drive_capacitance_ff": readout.drive_capacitance_ff,
        "critical_current_na": ic_a*1e9, "junction_area_um2": junction_area,
        "junction_linear_inductance_nh": PHI0/(2*math.pi*ic_a)*1e9,
        "effective_permittivity": eps_eff, "cpw_impedance_ohm": z0,
        "resonator_mode_capacitance_ff": mode_c/1e-15,
        "bare_half_wave_length_um": c/(2*fr*math.sqrt(eps_eff))*1e6,
        "external_q": q_external, "arm_reach_um": arm, "arm_width_um": width,
        "moat_gap_um": moat, "pad_scale": scale,
        "qubit_resonator_overlap_um": qr_overlap, "resonator_feed_overlap_um": rf_overlap,
        "drive_overlap_um": drive_overlap, "dispersive_estimate_valid": dispersive,
        "chi_mhz": chi_mhz, "purcell_t1_estimate_us": purcell_us,
    }


def wafer_sites(process: Process) -> list[dict]:
    radius = process.wafer_diameter_mm*500 - process.wafer_edge_exclusion_mm*1000
    if radius <= 0:
        raise ValueError("Wafer edge exclusion leaves no usable wafer area.")
    width, height = process.chip_width_um, process.chip_height_um
    px, py = width + process.dicing_street_um, height + process.dicing_street_um
    candidates = []
    for ox in [0, 0.5]:
        for oy in [0, 0.5]:
            sites = []
            for i in range(-math.ceil(radius/px)-1, math.ceil(radius/px)+2):
                for j in range(-math.ceil(radius/py)-1, math.ceil(radius/py)+2):
                    x, y = (i+ox)*px, (j+oy)*py
                    if math.hypot(abs(x)+width/2, abs(y)+height/2) <= radius:
                        sites.append({"x_um": x, "y_um": y})
            candidates.append(sites)
    result = max(candidates, key=lambda sites: (len(sites), -sum(p["x_um"]**2+p["y_um"]**2 for p in sites)))
    if not result:
        raise ValueError("No complete die fits inside the usable wafer area.")
    return [{"die": index+1, **site} for index, site in enumerate(sorted(result, key=lambda p:(p["y_um"], p["x_um"])))]


def _path(points, width, layer, extension=0):
    return gdstk.FlexPath(points, width, joins="natural", ends=(extension, extension), layer=layer).to_polygons()


def _serialize(polygons):
    return [[[round(float(x), 3), round(float(y), 3)] for x,y in polygon.points] for polygon in polygons]


def _csv(rows: list[dict]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def render_svg(geometry: dict, view: str = "chip") -> str:
    process, sizing = geometry["process"], geometry["sizing"]
    if view == "qubit":
        reach = max(sizing["arm_reach_um"]+150, 350)
        x0,x1,y0,y1 = -2300-reach, -2300+reach+250, -1300-reach, -1300+reach
    elif view == "junction":
        x = -2300
        y = -1300-sizing["arm_reach_um"]-sizing["moat_gap_um"]/2
        x0,x1,y0,y1 = x-15, x+15, y-15, y+15
    else:
        x0,x1 = -process["chip_width_um"]/2, process["chip_width_um"]/2
        y0,y1 = -process["chip_height_um"]/2, process["chip_height_um"]/2
    colors = {"ground":"#d8e0e9", "qubit":"#7c4ac0", "resonator":"#3478bd", "feedline":"#dc8b2c", "drive":"#239f7c", "al_qubit":"#db7093", "al_ground":"#db7093", "junction_intent":"#e64242"}
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0} {-y1} {x1-x0} {y1-y0}" role="img" aria-label="Preliminary chip layout in micrometres"><rect x="{x0}" y="{-y1}" width="{x1-x0}" height="{y1-y0}" fill="#f9fbfe"/><g transform="scale(1,-1)">']
    # Paint ground first so thin signal features remain visible at chip scale.
    for net in ["ground"]+[name for name in geometry["polygons_um"] if name!="ground"]:
        polygons=geometry["polygons_um"][net]
        for polygon in polygons:
            points = " ".join(f"{x},{y}" for x,y in polygon)
            pieces.append(f'<polygon points="{points}" fill="{colors[net]}" fill-rule="evenodd"/>')
    pieces.append('</g></svg>')
    return "".join(pieces)


def _dxf(geometry: dict) -> str:
    """Construction surfaces: ground blank minus explicit clearance polygons."""
    doc = ezdxf.new("R2010")
    doc.units = 13  # micrometres, DXF $INSUNITS
    model = doc.modelspace()
    process = geometry["process"]
    hw,hh = process["chip_width_um"]/2, process["chip_height_um"]/2
    layers = {"GROUND_BLANK": 8, "GROUND_SUBTRACT": 1, "NB_QUBIT": 6,
              "NB_RESONATOR": 5, "NB_FEEDLINE": 30, "NB_DRIVE": 3,
              "AL_QUBIT": 6, "AL_GROUND": 6, "JJ_INTENT_NOT_MASK": 1}
    for name,color in layers.items():
        doc.layers.new(name, dxfattribs={"color": color})
    model.add_lwpolyline([(-hw,-hh),(hw,-hh),(hw,hh),(-hw,hh)], close=True, dxfattribs={"layer":"GROUND_BLANK"})
    for polygon in geometry["ground_subtract_um"]:
        model.add_lwpolyline(polygon, close=True, dxfattribs={"layer":"GROUND_SUBTRACT"})
    names = {"qubit":"NB_QUBIT", "resonator":"NB_RESONATOR", "feedline":"NB_FEEDLINE", "drive":"NB_DRIVE", "al_qubit":"AL_QUBIT", "al_ground":"AL_GROUND", "junction_intent":"JJ_INTENT_NOT_MASK"}
    for net,name in names.items():
        for polygon in geometry["polygons_um"][net]:
            model.add_lwpolyline(polygon, close=True, dxfattribs={"layer":name})
    output = io.StringIO()
    doc.write(output)
    return output.getvalue()


def _spice(sizing: dict, readout: Readout) -> str:
    z0 = sizing["cpw_impedance_ohm"]
    return f"""* Linear small-signal equivalent; does not predict transmon anharmonicity.
* GHz, fF and nH parameters below are converted to SPICE units.
* Lossless half-wave line: physical end loading shifts the coupled modes.
VREAD read_src 0 AC 1
RIN read_src feed {z0:.9g}
ROUT feed 0 {z0:.9g}
CFEED feed r_out {sizing['resonator_feed_capacitance_ff']:.9g}f
TREAD r_out 0 r_in 0 Z0={z0:.9g} TD={1/(2*readout.frequency_ghz*1e9):.12g}
CQR r_in q {sizing['qubit_resonator_capacitance_ff']:.9g}f
CPAD q 0 {sizing['pad_capacitance_ff']:.9g}f
CJ q 0 {sizing['junction_capacitance_ff']:.9g}f
LJ q 0 {sizing['junction_linear_inductance_nh']:.9g}n
VDRIVE drive_src 0 AC 0
RDRIVE drive_src drive {z0:.9g}
CDRIVE drive q {readout.drive_capacitance_ff:.9g}f
.ac dec 200 1G 12G
.end
"""


EXPORT_NOTES = """PRELIMINARY DESIGN — EM EXTRACTION AND A PROCESS JUNCTION CELL ARE REQUIRED

Contents and coordinate conventions
----------------------------------
chip.gds: positive metal, hierarchical, units 1 um and grid 1 nm.
wafer.gds: repeated chip cell on a circular 2-inch wafer (diameter is configurable).
layout.svg: false-color preview of the same polygons, not a mask color map.
geometry.json: solver-neutral construction geometry, layers, materials, nets,
ports, circuit targets, analytical assumptions, and a content hash.
construction.dxf: micrometre construction surfaces, NOT an already-subtracted mask.
circuit.cir: linearized junction plus a lossless distributed resonator.
manifest.json, ports.csv, wafer_sites.csv, process.json: design provenance.
em_results_template.json: template for returning extracted values to the app.

GDS layers (check process.json for configured material-layer numbers)
------------------------------------------------------------------
Nb layer: ground, island, readout resonator, readout feedline, drive line.
Al layer: two landing leads overlapping Nb; no evaporation or contact recipe.
90: text annotations. 91: JJ target-area marker, NOT a fabrication junction.
92: port markers. 93: die outline. 94: wafer outline. 90–94 are not mask layers.
The Al leads intentionally stop at the JJ site. Insert your process-specific
junction cell, its overlaps/contacts, and relevant bridge/undercut layers.
The marker area comes from Ic/Jc; it is not a Dolan/Manhattan junction design.

HFSS / Q3D handoff
-----------------
Import construction.dxf in micrometres. Build ground as GROUND_BLANK minus
the UNION of GROUND_SUBTRACT polygons. Other NB_* and AL_* layers are conductor
surfaces. Extrude/assign film thicknesses from geometry.json if appropriate;
set silicon below z=0 and metal above z=0. Nb/Al names are identifiers, not
validated microwave material models. Set PEC, surface impedance or an appropriate
superconductor model explicitly. Do not import JJ_INTENT_NOT_MASK as metal.
Coordinates, reference impedance, and terminals are in ports.csv/geometry.json.
Assign Q3D conductor nets after uniting overlapping Nb/Al contact surfaces.
Define HFSS feed/drive excitations and a lumped JJ element between the two
jj_* terminals using the supplied LJ and CJ. Place terminal sheets/reference
grounds according to your solver setup; JSON port points alone are not ports.
This package has NOT been solved in HFSS/Q3D or run through a fabrication DRC.

Model limits
------------
CPW impedance and wavelength use a zero-thickness, dielectric half-space model;
substrate/film thicknesses are recorded for EM, not used in that approximation.
The half-wave length is the bare centerline length: bends, end capacitance,
coupler loading, kinetic inductance and mode hybridization shift its resonance.
Pad capacitance assumes uniform geometric scaling of an entered reference
cross. All three coupling capacitances use entered capacitance-per-overlap
calibrations at the specified gaps. Default calibrations are examples.
The target Csum budget includes pad, junction, drive and readout coupling;
effective EC still needs capacitance-matrix extraction for the coupled network.
g uses a weak-coupling transmon matrix-element approximation. chi and Purcell
estimates are only displayed when the entered parameters are dispersive;
the Purcell estimate is a two-level single-mode estimate, with no filter.
Airbridges, ground wirebonds, package modes, loss, process bias and flux lines
are not generated. Wafer packing ignores the wafer flat/notch and test structures.

Return extracted values
----------------------
Fill em_results_template.json with your own extraction values and source, keeping
geometry_sha256 unchanged. Load it in the app to compare targets and extracted
predictions. Uploaded numbers are user-supplied evidence, not an automatic solve.
"""


@dataclass
class LayoutBundle:
    geometry: dict
    manifest: dict
    chip_gds: bytes
    wafer_gds: bytes
    package_zip: bytes


def generate_layout(ej_ghz: float, ec_ghz: float, readout: Readout = Readout(),
                    process: Process = Process()) -> LayoutBundle:
    sizing = size_circuit(ej_ghz, ec_ghz, readout, process)
    layer = process.nb_layer
    qx,qy = -2300.0,-1300.0
    arm,width,moat = sizing["arm_reach_um"],sizing["arm_width_um"],sizing["moat_gap_um"]
    island = gdstk.boolean(
        gdstk.rectangle((qx-arm,qy-width/2),(qx+arm,qy+width/2)),
        gdstk.rectangle((qx-width/2,qy-arm),(qx+width/2,qy+arm)),
        "or", precision=0.001, layer=layer,
    )
    clearances = gdstk.offset(island, moat, join="miter", precision=0.001)
    metal = {"qubit": island}

    def add_cpw(name, points):
        metal[name] = _path(points, process.cpw_width_um, layer)
        clearances.extend(_path(points, process.cpw_width_um+2*process.cpw_gap_um, layer, process.cpw_gap_um))

    feed_y = 1800.0
    half_width,half_height = process.chip_width_um/2,process.chip_height_um/2
    inner_port_x = half_width-750
    add_cpw("feedline", [(-inner_port_x,feed_y),(inner_port_x,feed_y)])
    q_res_y = qy+width/2+process.qubit_resonator_gap_um+process.cpw_width_um/2
    start_x = qx+arm-sizing["qubit_resonator_overlap_um"]
    x_left = -1400.0
    res_feed_y = feed_y-process.cpw_width_um-process.resonator_feed_gap_um
    end_overlap = sizing["resonator_feed_overlap_um"]
    required_length = sizing["bare_half_wave_length_um"]
    fixed_length = (x_left-start_x)+(res_feed_y-q_res_y)+end_overlap
    route_width = None
    for lanes in [2,4,6]:
        span = (required_length-fixed_length)/lanes
        top_y = -500+(lanes-1)*350
        if 300 <= span <= 4000 and top_y+100 < res_feed_y:
            route_width = span
            break
    if route_width is None or start_x >= x_left or q_res_y >= -500:
        raise ValueError("The resonator target does not fit this routing template. Adjust targets or geometry calibration.")
    points = [(start_x,q_res_y),(x_left,q_res_y),(x_left,-500)]
    for lane in range(lanes):
        destination_x = x_left+route_width if lane%2==0 else x_left
        points.append((destination_x,-500+lane*350))
        if lane<lanes-1:
            points.append((destination_x,-500+(lane+1)*350))
    points.extend([(x_left,res_feed_y),(x_left+end_overlap,res_feed_y)])
    add_cpw("resonator",points)
    actual_length = sum(math.dist(a,b) for a,b in zip(points,points[1:]))
    if not math.isclose(actual_length,required_length,abs_tol=0.001):
        raise ValueError("The generated resonator length does not match its bare electrical length.")
    drive_y = qy-width/2-process.drive_gap_um-process.cpw_width_um/2
    drive_end = qx-arm+sizing["drive_overlap_um"]
    add_cpw("drive",[(-inner_port_x,drive_y),(drive_end,drive_y)])
    # Explicitly remove ground between the facing capacitor electrodes. Slot
    # envelopes alone would leave a screening ground strip at larger gaps.
    clearances.extend([
        gdstk.rectangle((start_x,qy+width/2),(qx+arm,q_res_y-process.cpw_width_um/2)),
        gdstk.rectangle((x_left,res_feed_y+process.cpw_width_um/2),
                        (x_left+end_overlap,feed_y-process.cpw_width_um/2)),
        gdstk.rectangle((qx-arm,drive_y+process.cpw_width_um/2),
                        (drive_end,qy-width/2)),
    ])

    def launch(net,y,side):
        # Pad, taper, and CPW are contiguous. Clearance is cut from the ground.
        a,b,d = side*(half_width-150),side*(half_width-450),side*inner_port_x
        pad = gdstk.rectangle((min(a,b),y-100),(max(a,b),y+100),layer=layer)
        taper = gdstk.Polygon([(b,y-100),(d,y-process.cpw_width_um/2),
                               (d,y+process.cpw_width_um/2),(b,y+100)],layer=layer)
        metal[net].extend([pad,taper])
        clearances.extend([
            gdstk.rectangle((min(a,b)-100,y-200),(max(a,b)+100,y+200)),
            gdstk.Polygon([(b,y-200),(d,y-process.cpw_width_um/2-process.cpw_gap_um),
                           (d,y+process.cpw_width_um/2+process.cpw_gap_um),(b,y+200)]),
        ])
        return (a+b)/2,y

    p_in = launch("feedline",feed_y,-1)
    p_out = launch("feedline",feed_y,1)
    p_drive = launch("drive",drive_y,-1)
    chip_rectangle = gdstk.rectangle((-half_width,-half_height),(half_width,half_height))
    metal["ground"] = gdstk.boolean(chip_rectangle,clearances,"not",precision=0.001,layer=layer)

    jj_y = qy-arm-moat/2
    jj_gap = 0.5  # Landing gap, NOT the fabricated junction barrier dimension.
    lead_width = max(1.0,process.minimum_nb_feature_um)
    metal["al_qubit"] = [gdstk.rectangle((qx-lead_width/2,jj_y+jj_gap/2),
                                         (qx+lead_width/2,qy-arm+1),layer=process.al_layer)]
    metal["al_ground"] = [gdstk.rectangle((qx-lead_width/2,qy-arm-moat-1),
                                          (qx+lead_width/2,jj_y-jj_gap/2),layer=process.al_layer)]
    jj_side = math.sqrt(sizing["junction_area_um2"])
    metal["junction_intent"] = [gdstk.rectangle((qx-jj_side/2,jj_y-jj_side/2),
                                                (qx+jj_side/2,jj_y+jj_side/2),layer=91)]

    # Check isolation of electrical nets, allowing intentional Nb/Al overlaps
    # only within a single net. The missing process-specific JJ stays open.
    nets = {"qubit":metal["qubit"]+metal["al_qubit"],
            "ground":metal["ground"]+metal["al_ground"],
            **{name:metal[name] for name in ["resonator","feedline","drive"]}}
    names = list(nets)
    for i,name in enumerate(names):
        for other in names[i+1:]:
            overlap = gdstk.boolean(nets[name],nets[other],"and",precision=0.001)
            if sum(p.area() for p in overlap)>1e-6:
                raise ValueError(f"Unintended metal short between {name} and {other}.")
    for name,polygons in metal.items():
        for polygon in polygons:
            box=polygon.bounding_box()
            if box[0][0]<-half_width-0.001 or box[1][0]>half_width+0.001 or box[0][1]<-half_height-0.001 or box[1][1]>half_height+0.001:
                raise ValueError(f"{name} is outside the chip boundary.")

    ports = [
        {"name":"readout_in","net":"feedline","x_um":p_in[0],"y_um":p_in[1],"reference_net":"ground","reference_ohm":sizing["cpw_impedance_ohm"]},
        {"name":"readout_out","net":"feedline","x_um":p_out[0],"y_um":p_out[1],"reference_net":"ground","reference_ohm":sizing["cpw_impedance_ohm"]},
        {"name":"drive_in","net":"drive","x_um":p_drive[0],"y_um":p_drive[1],"reference_net":"ground","reference_ohm":sizing["cpw_impedance_ohm"]},
        {"name":"jj_island_terminal","net":"qubit","x_um":qx,"y_um":jj_y+jj_gap/2,"reference_net":"ground","reference_ohm":0},
        {"name":"jj_ground_terminal","net":"ground","x_um":qx,"y_um":jj_y-jj_gap/2,"reference_net":"ground","reference_ohm":0},
    ]
    # Port coordinates must share the exported 1 nm grid, including the two
    # terminals at the edges of the narrow Al landing gap.
    for port in ports:
        port["x_um"]=round(port["x_um"],3)
        port["y_um"]=round(port["y_um"],3)
    geometry = {
        "schema_version":"transmon-layout/1.0", "units":"um", "status":"PRELIMINARY_NOT_EM_VERIFIED",
        "process":asdict(process), "readout_targets":asdict(readout), "sizing":sizing,
        "polygons_um":{name:_serialize(polygons) for name,polygons in metal.items()},
        "ground_subtract_um":_serialize(clearances), "ports":ports,
        "resonator_centerline_um":[[float(x),float(y)] for x,y in points],
        "resonator_drawn_length_um":actual_length,
        "z_stack_um":{"silicon":[-process.substrate_thickness_um,0],
                      "Nb":[0,process.nb_thickness_nm/1000],
                      "Al":[process.nb_thickness_nm/1000,(process.nb_thickness_nm+process.al_thickness_nm)/1000]},
    }
    canonical=json.dumps(geometry,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
    fingerprint=hashlib.sha256(canonical).hexdigest()
    geometry["geometry_sha256"]=fingerprint
    sites=wafer_sites(process)
    manifest={"schema_version":"transmon-layout/1.0","geometry_sha256":fingerprint,
              "status":"PRELIMINARY_NOT_EM_VERIFIED","calibration_source":process.calibration_source,
              "sizing":sizing,"process":asdict(process),"readout_targets":asdict(readout),
              "wafer_die_count":len(sites),
              "geometry_checks":{"no_unintended_net_overlap":True,"metal_inside_die":True,
                                 "bare_centerline_length_error_um":actual_length-required_length},
              "junction_process_cell_present":False,
              "ground_airbridges_and_package_not_included":True}
    library=gdstk.Library(name="TRANSMON_LAB",unit=1e-6,precision=1e-9)
    chip=library.new_cell("CHIP")
    for name,polygons in metal.items():
        child=library.new_cell(name.upper())
        child.add(*polygons)
        chip.add(gdstk.Reference(child))
    annotation=library.new_cell("ANNOTATIONS")
    annotation.add(gdstk.FlexPath([(-half_width,-half_height),(half_width,-half_height),
                                  (half_width,half_height),(-half_width,half_height),(-half_width,-half_height)],2,layer=93))
    for port in ports:
        annotation.add(gdstk.Label(port["name"],(port["x_um"],port["y_um"]+25),layer=90))
        annotation.add(gdstk.rectangle((port["x_um"]-5,port["y_um"]-5),(port["x_um"]+5,port["y_um"]+5),layer=92))
    chip.add(gdstk.Reference(annotation))
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        library.write_gds(root/"chip.gds",max_points=199)
        chip_data=(root/"chip.gds").read_bytes()
        wafer=library.new_cell("WAFER")
        for site in sites:
            wafer.add(gdstk.Reference(chip,origin=(site["x_um"],site["y_um"])))
        radius=process.wafer_diameter_mm*500
        wafer.add(gdstk.ellipse((0,0),radius,inner_radius=radius-10,tolerance=1,layer=94))
        library.write_gds(root/"wafer.gds",max_points=199)
        wafer_data=(root/"wafer.gds").read_bytes()
    output=io.BytesIO()
    em_template={"geometry_sha256":fingerprint,"source":"REPLACE with solver/run identifier",
                 "qubit_total_capacitance_ff":None,"readout_frequency_ghz":None,
                 "coupling_g_mhz":None,"loaded_linewidth_mhz":None}
    with ZipFile(output,"w",ZIP_DEFLATED) as archive:
        archive.writestr("chip.gds",chip_data)
        archive.writestr("wafer.gds",wafer_data)
        archive.writestr("layout.svg",render_svg(geometry))
        archive.writestr("construction.dxf",_dxf(geometry))
        archive.writestr("geometry.json",json.dumps(geometry,indent=2,allow_nan=False))
        archive.writestr("manifest.json",json.dumps(manifest,indent=2,allow_nan=False))
        archive.writestr("process.json",json.dumps(asdict(process),indent=2))
        archive.writestr("ports.csv",_csv(ports))
        archive.writestr("wafer_sites.csv",_csv(sites))
        archive.writestr("circuit.cir",_spice(sizing,readout))
        archive.writestr("em_results_template.json",json.dumps(em_template,indent=2))
        archive.writestr("README.txt",EXPORT_NOTES)
    return LayoutBundle(geometry,manifest,chip_data,wafer_data,output.getvalue())


def compare_extraction(bundle: LayoutBundle, results: dict) -> dict:
    if not isinstance(results,dict) or results.get("geometry_sha256")!=bundle.geometry["geometry_sha256"]:
        raise ValueError("The extraction file must have the geometry hash of the displayed layout.")
    if not isinstance(results.get("source"),str) or not results["source"].strip() or results["source"].startswith("REPLACE"):
        raise ValueError("Provide a solver/run identifier in the source field.")
    capacitance=results.get("qubit_total_capacitance_ff")
    if not isinstance(capacitance,(int,float)) or not math.isfinite(capacitance) or capacitance<=0:
        raise ValueError("Enter a positive extracted effective total qubit capacitance in fF.")
    predicted=observables(bundle.geometry["sizing"]["ej_ghz"],charging_energy_ghz(capacitance),charge_cutoff=45)
    compared={"source":results["source"],"extracted_f01_ghz":predicted.f01_ghz,
              "extracted_alpha_mhz":predicted.alpha_mhz,
              "frequency_error_mhz":(predicted.f01_ghz-bundle.geometry["sizing"]["target_f01_ghz"])*1000,
              "alpha_error_mhz":predicted.alpha_mhz-bundle.geometry["sizing"]["target_alpha_mhz"]}
    for name,target in [("readout_frequency_ghz",bundle.geometry["readout_targets"]["frequency_ghz"]),
                        ("coupling_g_mhz",bundle.geometry["readout_targets"]["coupling_g_mhz"]),
                        ("loaded_linewidth_mhz",bundle.geometry["readout_targets"]["loaded_linewidth_mhz"])]:
        value=results.get(name)
        if value is not None:
            if not isinstance(value,(float,int)) or not math.isfinite(value) or value<=0:
                raise ValueError(f"{name} must be a positive extracted value or null.")
            compared[name]=value
            compared[name+"_error"]=value-target
    return compared
