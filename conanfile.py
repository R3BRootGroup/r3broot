from conan import ConanFile
from conan.tools.cmake import CMakeToolchain


class CompressorRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"

    def requirements(self):
        self.requires("ms-gsl/4.0.0")
        self.requires("range-v3/0.12.0")
        self.requires("gtest/cci.20210126")
        self.requires("nlohmann_json/3.11.2")
        self.requires("yaml-cpp/0.8.0")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.user_presets_path = False
        tc.generate()
