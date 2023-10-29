# PICplusplus
This is a repo containing mostly C++ code for running Particle-In-Cell simulations of electrically-charged plasma distributions.

# How To Build
Note: use conan version < 2
1. conan install -pr <win_debug, win_release> -if <build_folder> --build=missing .
2. conan build -bf <build_folder> .
