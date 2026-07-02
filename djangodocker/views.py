from django.http import JsonResponse, HttpResponse
import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _find_executable() -> Path:
    candidates = [
        REPO_ROOT / "build" / "bin" / "PIC++Main",
        REPO_ROOT / "build" / "bin" / "Release" / "PIC++Main.exe",
        REPO_ROOT / "build" / "bin" / "PIC++Main.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "PIC++Main not found. Build the project first with ./scripts/build.sh"
    )


def _build_config(params: dict) -> dict:
    num_particles = int(params["numParticles"])
    drift_velocity = float(params["driftVelocity"])
    species_template = {
        "numParticles": num_particles,
        "spatialPerturbationMode": int(params["spatialPerturbationMode"]),
        "spatialPerturbationAmplitude": float(params["spatialPerturbationAmplitude"]),
        "thermalVelocity": float(params["thermalVelocity"]),
        "plasmaFrequency": float(params["plasmaFrequency"]),
        "chargeMassRatio": float(params["chargeMassRatio"]),
    }

    return {
        "species": [
            {"name": "Species1", "driftVelocity": drift_velocity, **species_template},
            {"name": "Species2", "driftVelocity": -drift_velocity, **species_template},
        ],
        "spatialLength": float(params["spatialLength"]),
        "numTimeSteps": int(params["timeSteps"]),
        "timeStepSize": float(params["timeStepSize"]),
        "numGrid": int(params["numGrid"]),
        "numSpecies": int(params["numSpecies"]),
    }


def run(request):
    logger.info("run view called.")

    source = request.GET if request.method == "GET" else request.POST
    required = [
        "spatialLength", "numParticles", "timeSteps", "timeStepSize", "numGrid",
        "spatialPerturbationMode", "driftVelocity", "numSpecies",
        "spatialPerturbationAmplitude", "thermalVelocity", "plasmaFrequency", "chargeMassRatio",
    ]

    missing = [key for key in required if source.get(key) is None]
    if missing:
        return HttpResponse(f"Missing parameters: {', '.join(missing)}", status=400)

    params = {key: source.get(key) for key in required}
    config = _build_config(params)

    input_path = None
    output_path = None

    try:
        executable = _find_executable()

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as input_file:
            json.dump(config, input_file)
            input_path = input_file.name

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as output_file:
            output_path = output_file.name

        result = subprocess.run(
            [str(executable), input_path, output_path],
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.error("PIC++Main failed: %s", result.stderr)
            return HttpResponse(
                f"Simulation failed:\n{result.stderr or result.stdout}",
                status=500,
            )

        with open(output_path, encoding="utf-8") as handle:
            data = json.load(handle)

        return JsonResponse(data)

    except FileNotFoundError as exc:
        return HttpResponse(str(exc), status=500)
    except json.JSONDecodeError as exc:
        return HttpResponse(f"Invalid simulation output: {exc}", status=500)
    except Exception as exc:
        logger.exception("Error running simulation")
        return HttpResponse(f"Error executing C++ program: {exc}", status=500)
    finally:
        for path in (input_path, output_path):
            if path and os.path.exists(path):
                os.unlink(path)
