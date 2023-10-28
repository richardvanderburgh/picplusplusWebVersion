from conans import ConanFile, CMake
from conans.tools import download, unzip, check_md5, check_sha1, check_sha256
import os
import shutil

class pic_plus_plus(ConanFile):
    name = "PIC++"
    version = "0.0.1"
    generators = "cmake", "cmake_find_package"

    def requirements(self):
        self.requires("gtest/1.13.0")
        self.requires("nlohmann_json/3.11.2")

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        # here you can run CTest, launch your binaries, etc
        cmake.test(target="RUN_TESTS", output_on_failure=True)