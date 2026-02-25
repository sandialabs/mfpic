# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack.package import *

class Mfpic(CMakePackage):
    """Multi-Fidelity Particle-In-Cell."""

    homepage = "https://github.com/sandialabs/mfpic"
    git = "https://github.com/sandialabs/mfpic.git"

    maintainers(
        "reuterb"
        "logantm2",
        "caleblogemann",
        "gyalla",
    )

    license("BSD-3-Clause", checked_by="logantm2")

    version("main", branch="main")

    generator("ninja")

    depends_on("cxx", type="build")
    depends_on("c",   type="build")
    depends_on("cmake@3.10:", type="build")
    depends_on("googletest", type="test")
    depends_on("mfem~zlib~mpi~metis@=4.8.0", when="build_type=Release")
    depends_on("mfem~zlib~mpi~metis@=4.8.0 +debug", when="build_type=Debug")
    depends_on("ninja", type="build")
    depends_on("python@3:", type="test")
    depends_on("py-numpy", type="test")
    depends_on("py-scipy", type="test")
    depends_on("vtk~mpi~opengl2+python", type="test")
    depends_on("py-matplotlib", type="test")
    depends_on("yaml-cpp")
    depends_on("hdf5+hl", when="build_type=Release")
    depends_on("hdf5+hl build_type=Debug", when="build_type=Debug")
    depends_on("gdb", when="build_type=Debug")

    @property
    def keep_werror(self):
        return "all"

    requires(
        "%gcc@11:",
        "%clang@14:",
        policy="one_of",
        msg="builds only with a GCC or Clang that supports C++ 23",
    )

    requires(
        "^mfem~static",
        when="%llvm",
        msg="Clang compiler requires PIC mfem",
    )
