"""
URL configuration for PICplusplusServices project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from .views import run, run_start, run_status, run_result, demo_list, demo_detail

from django.http import HttpResponse
from django.views.generic import TemplateView

def home(request):
    return HttpResponse("Welcome to my Django app!")

urlpatterns = [
        path('', TemplateView.as_view(template_name="index.html")),
        path('run/', run, name='run'),
        path('run/start/', run_start, name='run_start'),
        path('run/status/<str:job_id>/', run_status, name='run_status'),
        path('run/result/<str:job_id>/', run_result, name='run_result'),
        path('demos/', demo_list, name='demo_list'),
        path('demos/<str:demo_id>/', demo_detail, name='demo_detail'),
]
