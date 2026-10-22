from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_explore_and_design_user_workflow():
    app = AppTest.from_file(str(Path(__file__).parents[1] / "app.py"), default_timeout=30).run()
    assert not app.exception
    initial_frequency = next(x.value for x in app.metric if x.label == "Qubit frequency · f01")
    app.slider(key="capacitance_ff").set_value(120).run()
    assert not app.exception
    changed_frequency = next(x.value for x in app.metric if x.label == "Qubit frequency · f01")
    assert float(changed_frequency.split()[0]) < float(initial_frequency.split()[0])

    app.number_input(key="target_f01").set_value(6.2)
    app.number_input(key="target_alpha").set_value(-200)
    app.button[0].click().run()
    assert not app.exception
    assert len(app.success) == 1
    assert "6.2000 GHz" in app.success[0].value
    assert "-200.00 MHz" in app.success[0].value

    app.slider(key="ej_ghz").set_value(1).run()
    assert not app.exception
    assert any("below 20" in warning.value for warning in app.warning)
    # An unrelated exploration change must not relabel the saved design target.
    assert "6.2000 GHz" in app.success[0].value


def test_target_to_gds_user_workflow():
    app=AppTest.from_file(str(Path(__file__).parents[1]/"app.py"),default_timeout=30).run()
    app.button[0].click().run()
    assert not app.exception
    assert app.selectbox(key="layout_source").value=="Saved target design"
    next(button for button in app.button if button.label=="Generate chip layout").click().run()
    assert not app.exception
    assert not app.error
    bundle=app.session_state["layout_bundle"]
    assert abs(bundle.geometry["sizing"]["target_f01_ghz"]-5)<1e-6
    assert len(bundle.chip_gds)>1000
    assert bundle.manifest["wafer_die_count"]==12
    app.selectbox(key="layout_view").set_value("junction").run()
    assert not app.exception
