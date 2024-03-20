from conan import ConanFile
from conan.tools.cmake import cmake_layout, CMake, CMakeToolchain
import os
import shutil

class pic_plus_plus(ConanFile):
    name = "PIC++"
    version = "0.0.1"
    generators = "CMakeDeps"
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        self.requires("gtest/1.13.0")
        self.requires("nlohmann_json/3.11.2")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = "ON"
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        # here you can run CTest, launch your binaries, etc
        cmake.test()

    
