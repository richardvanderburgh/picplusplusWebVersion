from django.http import JsonResponse, HttpResponse
import json
import logging
import os
import subprocess
import tempfile
import threading
import uuid
from pathlib import Path

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEMO_DIR = REPO_ROOT / "inputFiles" / "demo"

JOBS = {}
JOBS_LOCK = threading.Lock()


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


def _parse_params(source) -> dict:
    required = [
        "spatialLength", "numParticles", "timeSteps", "timeStepSize", "numGrid",
        "spatialPerturbationMode", "driftVelocity", "numSpecies",
        "spatialPerturbationAmplitude", "thermalVelocity", "plasmaFrequency",
        "chargeMassRatio",
    ]
    missing = [key for key in required if source.get(key) is None]
    if missing:
        raise ValueError(f"Missing parameters: {', '.join(missing)}")
    params = {key: source.get(key) for key in required}
    params["spatialPerturbationWaveform"] = source.get("spatialPerturbationWaveform", "cos")
    params["framePeriod"] = source.get("framePeriod", "1")

    num_species = int(params["numSpecies"])
    if num_species not in (1, 2):
        raise ValueError("numSpecies must be 1 or 2")

    return params


def _build_config(params: dict) -> dict:
    num_particles = int(params["numParticles"])
    num_species = int(params["numSpecies"])
    drift_velocity = float(params["driftVelocity"])

    species_template = {
        "numParticles": num_particles,
        "spatialPerturbationMode": int(params["spatialPerturbationMode"]),
        "spatialPerturbationAmplitude": float(params["spatialPerturbationAmplitude"]),
        "spatialPerturbationWaveform": str(params["spatialPerturbationWaveform"]),
        "thermalVelocity": float(params["thermalVelocity"]),
        "plasmaFrequency": float(params["plasmaFrequency"]),
        "chargeMassRatio": float(params["chargeMassRatio"]),
    }

    if num_species == 1:
        species = [
            {"name": "Species1", "driftVelocity": drift_velocity, **species_template},
        ]
    else:
        species = [
            {"name": "Species1", "driftVelocity": drift_velocity, **species_template},
            {"name": "Species2", "driftVelocity": -drift_velocity, **species_template},
        ]

    frame_period = int(params["framePeriod"])
    if frame_period < 1:
        frame_period = 1

    return {
        "species": species,
        "spatialLength": float(params["spatialLength"]),
        "numTimeSteps": int(params["timeSteps"]),
        "timeStepSize": float(params["timeStepSize"]),
        "numGrid": int(params["numGrid"]),
        "numSpecies": len(species),
        "framePeriod": frame_period,
    }


def _set_job(job_id: str, **fields) -> None:
    with JOBS_LOCK:
        JOBS.setdefault(job_id, {})
        JOBS[job_id].update(fields)


def _run_job(job_id: str, config: dict, input_path: str, output_path: str) -> None:
    try:
        executable = _find_executable()
        _set_job(job_id, status="running", percent=0, message="Starting PIC++Main")

        env = os.environ.copy()
        env["PICPP_PROGRESS"] = "1"

        with open(input_path, "w", encoding="utf-8") as handle:
            json.dump(config, handle)

        process = subprocess.Popen(
            [str(executable), input_path, output_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )

        assert process.stderr is not None
        stderr_tail = []
        for line in process.stderr:
            line = line.strip()
            stderr_tail.append(line)
            if line.startswith("PICPP_PROGRESS "):
                try:
                    percent = int(line.split()[1])
                    percent = max(0, min(100, percent))
                    _set_job(job_id, percent=percent, message=f"Time loop {percent}%")
                except (IndexError, ValueError):
                    logger.warning("Unparseable progress line: %s", line)

        returncode = process.wait()
        if returncode != 0:
            _set_job(
                job_id,
                status="error",
                percent=0,
                message="Simulation failed",
                error="\n".join(stderr_tail[-20:]) or f"exit code {returncode}",
            )
            return

        with open(output_path, encoding="utf-8") as handle:
            result = json.load(handle)

        _set_job(job_id, status="complete", percent=100, message="Done", result=result)
    except Exception as exc:
        logger.exception("Job %s failed", job_id)
        _set_job(job_id, status="error", message="Simulation failed", error=str(exc))
    finally:
        for path in (input_path, output_path):
            try:
                if path and os.path.exists(path):
                    os.unlink(path)
            except OSError:
                pass


def run_start(request):
    source = request.GET if request.method == "GET" else request.POST
    try:
        params = _parse_params(source)
        config = _build_config(params)
    except ValueError as exc:
        return HttpResponse(str(exc), status=400)

    job_id = str(uuid.uuid4())

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as input_file:
        input_path = input_file.name
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as output_file:
        output_path = output_file.name

    _set_job(
        job_id,
        status="queued",
        percent=0,
        message="Queued",
        error=None,
        result=None,
    )

    thread = threading.Thread(
        target=_run_job,
        args=(job_id, config, input_path, output_path),
        daemon=True,
    )
    thread.start()

    return JsonResponse({"jobId": job_id})


def run_status(request, job_id):
    with JOBS_LOCK:
        job = JOBS.get(job_id)
    if job is None:
        return JsonResponse({"error": "Unknown job"}, status=404)
    return JsonResponse(
        {
            "status": job.get("status", "unknown"),
            "percent": job.get("percent", 0),
            "message": job.get("message", ""),
            "error": job.get("error"),
        }
    )


def run_result(request, job_id):
    with JOBS_LOCK:
        job = JOBS.get(job_id)
        if job is None:
            return JsonResponse({"error": "Unknown job"}, status=404)
        if job.get("status") == "error":
            return JsonResponse({"error": job.get("error", "Simulation failed")}, status=500)
        if job.get("status") != "complete" or job.get("result") is None:
            return JsonResponse({"error": "Job not complete"}, status=409)
        result = job["result"]
        job.pop("result", None)

    return JsonResponse(result)


def run(request):
    """Synchronous run (legacy). Prefer /run/start/ + polling for progress."""
    source = request.GET if request.method == "GET" else request.POST
    try:
        params = _parse_params(source)
        config = _build_config(params)
    except ValueError as exc:
        return HttpResponse(str(exc), status=400)

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


def _load_demo_manifest() -> list[dict]:
    manifest_path = DEMO_DIR / "manifest.json"
    with open(manifest_path, encoding="utf-8") as handle:
        return json.load(handle)


def _load_demo_config(demo_id: str) -> dict:
    for entry in _load_demo_manifest():
        if entry["id"] == demo_id:
            demo_path = DEMO_DIR / entry["file"]
            with open(demo_path, encoding="utf-8") as handle:
                return json.load(handle)
    raise FileNotFoundError(f"Unknown demo: {demo_id}")


def config_to_form_params(config: dict) -> dict:
    if not config.get("species"):
        raise ValueError("Demo config has no species")
    species0 = config["species"][0]
    frame_period = config.get("framePeriod", 5)
    if frame_period < 1:
        frame_period = 5

    return {
        "spatialLength": str(config["spatialLength"]),
        "numParticles": str(species0["numParticles"]),
        "timeSteps": str(config["numTimeSteps"]),
        "timeStepSize": str(config["timeStepSize"]),
        "numGrid": str(config["numGrid"]),
        "spatialPerturbationMode": str(species0.get("spatialPerturbationMode", 0)),
        "driftVelocity": str(species0.get("driftVelocity", 0)),
        "numSpecies": str(config["numSpecies"]),
        "spatialPerturbationAmplitude": str(species0.get("spatialPerturbationAmplitude", 0)),
        "thermalVelocity": str(species0.get("thermalVelocity", 0)),
        "plasmaFrequency": str(species0.get("plasmaFrequency", 1)),
        "chargeMassRatio": str(species0.get("chargeMassRatio", -1)),
        "spatialPerturbationWaveform": species0.get("spatialPerturbationWaveform", "cos"),
        "framePeriod": str(frame_period),
    }


def demo_list(request):
    try:
        demos = _load_demo_manifest()
    except (OSError, json.JSONDecodeError) as exc:
        return JsonResponse({"error": str(exc)}, status=500)
    return JsonResponse({"demos": demos})


def demo_detail(request, demo_id):
    try:
        config = _load_demo_config(demo_id)
        params = config_to_form_params(config)
        entry = next(item for item in _load_demo_manifest() if item["id"] == demo_id)
    except StopIteration:
        return JsonResponse({"error": "Unknown demo"}, status=404)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        return JsonResponse({"error": str(exc)}, status=500)

    return JsonResponse(
        {
            "id": entry["id"],
            "title": entry["title"],
            "category": entry.get("category", ""),
            "description": entry.get("description", ""),
            "lookFor": entry.get("lookFor", []),
            "params": params,
            "config": config,
        }
    )
