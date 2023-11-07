# PICplusplus
This is a repo containing a C++ executable and HTML UI for running Particle-In-Cell simulations of electrically-charged plasma distributions.

# How To Build
Note: use conan version >= 2
1. conan install -pr:h=<win_debug, win_release> -pr:b=<win_debug, win_release> -of <build_folder> --build=missing .
2. conan build -of <build_folder> .

# How to Run the UI in a Browser
Note: requires Python and Django 
1. python manage.py runserver
2. Open the development server