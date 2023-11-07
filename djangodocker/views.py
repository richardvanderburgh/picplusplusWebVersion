from django.http import JsonResponse
import subprocess
from django.http import HttpResponse
import json

import logging

logger = logging.getLogger(__name__)

def run(request):
    logger.info('run view called.')

    if request.method == 'GET':
        spatialLength =  request.GET.get('spatialLength')
        numParticles = request.GET.get('numParticles')
        timeSteps = request.GET.get('timeSteps')       
        timeStepSize = request.GET.get('timeStepSize')
        numGrid = request.GET.get('numGrid')
        spatialPerturbationMode = request.GET.get('spatialPerturbationMode')
        driftVelocity = request.GET.get('driftVelocity')
        numSpecies = request.GET.get('numSpecies')
        spatialPerturbationAmplitude = request.GET.get('spatialPerturbationAmplitude')
        thermalVelocity = request.GET.get('thermalVelocity')
        plasmaFrequency = request.GET.get('plasmaFrequency')
        chargeMassRatio= request.GET.get('chargeMassRatio')

    else:
        spatialLength =  request.POST.get('spatialLength')
        numParticles = request.POST.get('numParticles')
        timeSteps = request.POST.get('timeSteps')
        timeStepSize = request.POST.get('timeStepSize')
        numGrid = request.POST.get('numGrid')
        spatialPerturbationMode = request.POST.get('spatialPerturbationMode')
        driftVelocity = request.POST.get('driftVelocity')
        numSpecies = request.POST.get('numSpecies')
        spatialPerturbationAmplitude = request.POST.get('spatialPerturbationAmplitude')
        thermalVelocity = request.POST.get('thermalVelocity')
        plasmaFrequency = request.POST.get('plasmaFrequency')
        chargeMassRatio= request.POST.get('chargeMassRatio')

    # Run the C++ executable
    try:
        executablePath = "build/bin/Release/PIC++Main.exe"

        result = subprocess.run([executablePath, spatialLength, numParticles, timeSteps, timeStepSize, numGrid, spatialPerturbationMode, driftVelocity, numSpecies, spatialPerturbationAmplitude, thermalVelocity, plasmaFrequency, chargeMassRatio], capture_output=True, text=True)
        output = result.stdout.strip()

        return JsonResponse(output, safe=False)
        #return HttpResponse("C++ program executed successfully.")
    except Exception as e:
        return HttpResponse(f"Error executing C++ program: {str(e)}, {executablePath}")
