import os
import subprocess

from setuptools import setup, Distribution
from setuptools.command.build_py import build_py


class BuildWithExecutable(build_py):
    """Normal Python build, then compile the C solver for this machine."""

    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return  # developers compile by hand: cd wind_ae/src && make compile
        self._compile_executable()

    def _compile_executable(self):
        here = os.path.abspath(os.path.dirname(__file__))
        src_dir = os.path.join(here, "wind_ae", "src")
        out_dir = os.path.abspath(os.path.join(self.build_lib, "wind_ae", "bin"))
        os.makedirs(out_dir, exist_ok=True)
        dest = os.path.join(out_dir, "relaxed_ae")
        try:
            subprocess.run(["make", "compile", "BIN=" + dest], cwd=src_dir, check=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as err:
            raise RuntimeError(
                "Wind-AE could not compile its C solver (relaxed_ae). "
                "You need `make` and a C compiler. Original error: " + str(err)
            )
        if not os.path.exists(dest):
            raise RuntimeError("make finished but relaxed_ae was not produced.")
        os.chmod(dest, 0o755)


class PlatformDistribution(Distribution):
    """Tag wheels built on the user's machine as platform-specific."""

    def has_ext_modules(self):
        return True


setup(cmdclass={"build_py": BuildWithExecutable}, distclass=PlatformDistribution)
