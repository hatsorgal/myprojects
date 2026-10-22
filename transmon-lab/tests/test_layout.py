from dataclasses import asdict, replace
import hashlib
import io
import json
import math
from zipfile import ZipFile

import ezdxf
import gdstk
import numpy as np
import pytest
from scipy.constants import elementary_charge, h

from chip_layout import (
    Process, Readout, compare_extraction, cpw_parameters, generate_layout,
    load_process, render_svg, size_circuit, wafer_sites,
)
from transmon import design_transmon


@pytest.fixture(scope="module")
def bundle():
    target=design_transmon(5,-250)
    return generate_layout(target.ej_ghz,target.ec_ghz)


def test_cpw_and_circuit_units(bundle):
    eps,z0=cpw_parameters(Process())
    assert eps==pytest.approx(6.225)
    assert 45<z0<55
    sizing=bundle.geometry["sizing"]
    assert 8000<sizing["bare_half_wave_length_um"]<9000
    expected_ic=4*math.pi*elementary_charge*sizing["ej_ghz"]*1e9
    assert sizing["critical_current_na"]==pytest.approx(expected_ic*1e9)
    assert sum(sizing[name] for name in ["pad_capacitance_ff","junction_capacitance_ff",
                                        "qubit_resonator_capacitance_ff","drive_capacitance_ff"])==pytest.approx(sizing["total_capacitance_ff"])
    assert sizing["chi_mhz"]<0
    assert sizing["purcell_t1_estimate_us"]>0


def test_gds_roundtrip_units_layers_and_resonator_length(bundle,tmp_path):
    path=tmp_path/"chip.gds"
    path.write_bytes(bundle.chip_gds)
    library=gdstk.read_gds(path)
    assert library.unit==pytest.approx(1e-6)
    assert library.precision==pytest.approx(1e-9)
    cells={cell.name:cell for cell in library.cells}
    assert {cell.name for cell in library.top_level()}=={"CHIP"}
    assert {1,2,91,92,93} <= {p.layer for cell in library.cells for p in cell.get_polygons()}
    res_area=sum(p.area() for p in cells["RESONATOR"].get_polygons())
    expected=bundle.geometry["resonator_drawn_length_um"]*bundle.geometry["process"]["cpw_width_um"]
    # Export grid rounding is 1 nm; allow area error from rounded vertices.
    assert res_area==pytest.approx(expected,abs=0.1)
    assert bundle.manifest["junction_process_cell_present"] is False
    assert "AL_QUBIT" in cells and "AL_GROUND" in cells


def test_exported_nets_are_disjoint_and_ports_lie_on_their_nets(bundle):
    shapes=bundle.geometry["polygons_um"]
    nets={"qubit":shapes["qubit"]+shapes["al_qubit"],
          "ground":shapes["ground"]+shapes["al_ground"],
          **{name:shapes[name] for name in ["resonator","feedline","drive"]}}
    names=list(nets)
    for i,name in enumerate(names):
        for other in names[i+1:]:
            intersection=gdstk.boolean(nets[name],nets[other],"and",precision=0.001)
            assert sum(p.area() for p in intersection)==0,(name,other)
    for port in bundle.geometry["ports"]:
        assert gdstk.inside([(port["x_um"],port["y_um"])],nets[port["net"]])[0]


def test_wafer_all_die_corners_fit_and_streets_are_clear(bundle,tmp_path):
    process=Process()
    sites=wafer_sites(process)
    assert len(sites)==12
    radius=process.wafer_diameter_mm*500-process.wafer_edge_exclusion_mm*1000
    for i,site in enumerate(sites):
        assert math.hypot(abs(site["x_um"])+5000,abs(site["y_um"])+5000)<=radius
        for other in sites[i+1:]:
            assert abs(site["x_um"]-other["x_um"])>=10500 or abs(site["y_um"]-other["y_um"])>=10500
    path=tmp_path/"wafer.gds"
    path.write_bytes(bundle.wafer_gds)
    library=gdstk.read_gds(path)
    top=library.top_level()
    assert len(top)==1 and top[0].name=="WAFER"
    assert len(top[0].references)==len(sites)


def test_solver_package_integrity_dxf_units_and_geometry_hash(bundle):
    with ZipFile(io.BytesIO(bundle.package_zip)) as archive:
        assert archive.testzip() is None
        assert {"chip.gds","wafer.gds","geometry.json","construction.dxf","circuit.cir",
                "ports.csv","process.json","em_results_template.json","README.txt"}<=set(archive.namelist())
        geometry=json.loads(archive.read("geometry.json"))
        fingerprint=geometry.pop("geometry_sha256")
        canonical=json.dumps(geometry,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
        assert hashlib.sha256(canonical).hexdigest()==fingerprint
        drawing=ezdxf.read(io.StringIO(archive.read("construction.dxf").decode()))
        assert drawing.units==13
        assert not drawing.audit().has_errors
        names={entity.dxf.layer for entity in drawing.modelspace()}
        assert {"GROUND_BLANK","GROUND_SUBTRACT","NB_QUBIT","NB_RESONATOR","AL_QUBIT"}<=names
        assert "Linear small-signal" in archive.read("circuit.cir").decode()


def test_extraction_rejects_a_different_geometry(bundle):
    with pytest.raises(ValueError,match="geometry hash"):
        compare_extraction(bundle,{"geometry_sha256":"stale-layout"})


def test_extraction_predicts_effect_of_changed_capacitance(bundle):
    data={"geometry_sha256":bundle.geometry["geometry_sha256"],"source":"Q3D test fixture",
          "qubit_total_capacitance_ff":bundle.geometry["sizing"]["total_capacitance_ff"],
          "readout_frequency_ghz":6.9}
    compared=compare_extraction(bundle,data)
    assert compared["extracted_f01_ghz"]==pytest.approx(5,abs=1e-6)
    assert compared["readout_frequency_ghz_error"]==pytest.approx(-0.1)
    data["qubit_total_capacitance_ff"]*=1.1
    assert compare_extraction(bundle,data)["extracted_f01_ghz"]<5


def test_coupler_window_has_no_screening_ground_at_a_larger_gap():
    process=replace(Process(),qubit_resonator_gap_um=30,resonator_feed_gap_um=30,drive_gap_um=30)
    target=design_transmon(5,-250)
    bundle=generate_layout(target.ej_ghz,target.ec_ghz,process=process)
    sizing=bundle.geometry["sizing"]
    x=-2300+sizing["arm_reach_um"]-sizing["qubit_resonator_overlap_um"]/2
    y=-1300+sizing["arm_width_um"]/2+15
    assert not gdstk.inside([(x,y)],bundle.geometry["polygons_um"]["ground"])[0]


def test_invalid_coupling_budget_or_intrinsic_linewidth_rejected():
    with pytest.raises(ValueError,match="intrinsic linewidth"):
        size_circuit(20,0.24,Readout(loaded_linewidth_mhz=0.01,intrinsic_q=1000),Process())
    with pytest.raises(ValueError,match="does not fit"):
        generate_layout(20,0.24,Readout(coupling_g_mhz=150))


def test_non_dispersive_estimates_are_withheld():
    fq=5.0
    target=design_transmon(fq,-250)
    sizing=size_circuit(target.ej_ghz,target.ec_ghz,Readout(frequency_ghz=fq),Process())
    assert sizing["dispersive_estimate_valid"] is False
    assert sizing["chi_mhz"] is None and sizing["purcell_t1_estimate_us"] is None


@pytest.mark.parametrize("changes",[{"nb_layer":2},{"substrate_eps_r":math.nan},{"unrecognized":4},
                                    {"wafer_diameter_mm":100000},{"substrate_thickness_um":1},
                                    {"junction_jc_a_cm2":True}])
def test_invalid_process_values_rejected(changes):
    with pytest.raises(ValueError):
        load_process({**asdict(Process()),**changes})


def test_all_svg_views_are_self_contained_xml(bundle):
    from xml.etree import ElementTree
    for view in ["chip","qubit","junction"]:
        root=ElementTree.fromstring(render_svg(bundle.geometry,view))
        assert root.tag.endswith("svg")
        assert root.findall(".//{http://www.w3.org/2000/svg}polygon")
