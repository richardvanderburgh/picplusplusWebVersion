#!/usr/bin/env python3
"""Smoke tests for Django config builders (no server required)."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from djangodocker.views import (  # noqa: E402
    _build_config,
    _validate_config,
    config_to_form_params,
)


def test_two_species_opposite_drift():
    params = {
        "spatialLength": "6.28",
        "numParticles": "100",
        "timeSteps": "10",
        "timeStepSize": "0.1",
        "numGrid": "32",
        "spatialPerturbationMode": "1",
        "driftVelocity": "1.5",
        "numSpecies": "2",
        "spatialPerturbationAmplitude": "0.01",
        "thermalVelocity": "0",
        "plasmaFrequency": "1",
        "chargeMassRatio": "-1",
        "spatialPerturbationWaveform": "cos",
        "framePeriod": "1",
    }
    config = _build_config(params)
    assert config["numSpecies"] == 2
    assert config["species"][0]["driftVelocity"] == 1.5
    assert config["species"][1]["driftVelocity"] == -1.5


def test_full_config_preserves_asymmetric_species():
    config = {
        "species": [
            {
                "name": "A",
                "numParticles": 500,
                "driftVelocity": 0.0,
                "thermalVelocity": 0.5,
                "spatialPerturbationMode": 1,
                "spatialPerturbationAmplitude": 0.02,
                "plasmaFrequency": 1,
                "chargeMassRatio": -1,
            },
            {
                "name": "B",
                "numParticles": 300,
                "driftVelocity": 0.1,
                "thermalVelocity": 0.4,
                "spatialPerturbationMode": 2,
                "spatialPerturbationAmplitude": 0.01,
                "plasmaFrequency": 2,
                "chargeMassRatio": -2,
            },
        ],
        "spatialLength": 6.28318530717958,
        "numTimeSteps": 300,
        "timeStepSize": 0.1,
        "numGrid": 256,
        "numSpecies": 2,
    }
    validated = _validate_config(config)
    assert validated["species"][0]["thermalVelocity"] == 0.5
    assert validated["species"][1]["chargeMassRatio"] == -2
    assert validated["species"][1]["numParticles"] == 300

    form = config_to_form_params(validated)
    assert form["numSpecies"] == "2"
    assert form["numGrid"] == "256"
    assert form["thermalVelocity"] == "0.5"


def test_rejects_non_power_of_two_grid():
    try:
        _validate_config(
            {
                "species": [{"name": "A", "numParticles": 10}],
                "spatialLength": 1.0,
                "numTimeSteps": 1,
                "timeStepSize": 0.1,
                "numGrid": 30,
                "numSpecies": 1,
            }
        )
    except ValueError as exc:
        assert "power of two" in str(exc)
    else:
        raise AssertionError("expected ValueError for numGrid=30")


def test_example_input_json_round_trip():
    import json

    path = ROOT / "inputFiles" / "exampleInput.json"
    config = json.loads(path.read_text())
    validated = _validate_config(config)
    assert validated["species"][0]["thermalVelocity"] == 0.5
    assert validated["species"][1]["thermalVelocity"] == 0.4


if __name__ == "__main__":
    test_two_species_opposite_drift()
    test_full_config_preserves_asymmetric_species()
    test_rejects_non_power_of_two_grid()
    test_example_input_json_round_trip()
    print("OK: django config smoke tests passed")
