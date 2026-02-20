# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install mfpic
#
# You can edit this file again by typing:
#
#     spack edit mfpic
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack.package import *


class Mfpic(CMakePackage):
    """FIXME: Put a proper description of your package here."""

    homepage = "https://github.com/sandialabs/mfpic"
    git = "git@github.com:sandialabs/mfpic.git"

    # FIXME: Add a list of GitHub accounts to
    # notify when the package is updated.
    # maintainers("github_user1", "github_user2")

    # FIXME: Add the SPDX identifier of the project's license below.
    # See https://spdx.org/licenses/ for a list. Upon manually verifying
    # the license, set checked_by to your Github username.
    license("UNKNOWN", checked_by="github_user1")

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
    depends_on("yaml-cpp")
    depends_on("hdf5", when="build_type=Release")
    depends_on("hdf5 build_type=Debug", when="build_type=Debug")
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
