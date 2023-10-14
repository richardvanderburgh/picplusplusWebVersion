from django.http import JsonResponse
import subprocess
from django.http import HttpResponse
import json

import logging

logger = logging.getLogger(__name__)

def run_simulation(request):
    logger.info('run_simulation view called.')

    if request.method == 'GET':
        param1 = request.GET.get('param1')
        param2 = request.GET.get('param2')       
        param3 = request.GET.get('param3')
        param4 = request.GET.get('param4')
        param5 = request.GET.get('param5')
        param6 = request.GET.get('param6')
        param7 = request.GET.get('param7')
        param8 = request.GET.get('param8')


    else:
        param1 = request.POST.get('param1')
        param2 = request.POST.get('param2')
        param3 = request.POST.get('param3')
        param4 = request.POST.get('param4')
        param5 = request.POST.get('param5')
        param6 = request.POST.get('param6')
        param7 = request.POST.get('param7')
        param8 = request.POST.get('param8')

    # Run the C++ executable
    try:
        executablePath = "/PIC++Main"
        result = subprocess.run([executablePath, param1, param2, param3, param4, param5, param6, param7, param8], capture_output=True, text=True)
        output = result.stdout.strip()

        return JsonResponse(output, safe=False)
        #return HttpResponse("C++ program executed successfully.")
    except Exception as e:
        return HttpResponse(f"Error executing C++ program: {str(e)}")
