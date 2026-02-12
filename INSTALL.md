# For developers

1. Install [Spack prerequisites](https://spack.readthedocs.io/en/latest/installing_prerequisites.html).
   On Ubuntu/Debian:
   ```bash
   sudo apt update
   sudo apt install file bzip2 ca-certificates g++ gcc gfortran git gzip lsb-release patch python3 tar unzip xz-utils zstd
   ```

1. Install [Spack](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html). Usually:
   ```bash
   git clone --depth=2 https://github.com/spack/spack.git ~/spack
   ```
   Spack stores a bunch of user-local cache and config stuff in `~/.spack`.
   To override this, you can add
   ```bash
   export SPACK_DISABLE_LOCAL_CONFIG=true
   export SPACK_USER_CACHE_PATH=~/spack_cache
   ```
   to your `.bashrc`. Then do an initial Spack setup with
   ```bash
   . ~/spack/share/spack/setup-env.sh # Recommend adding this line to your .bashrc
   spack bootstrap now
   spack compiler find
   ```

1. Create and enter a [Spack environment](https://spack-tutorial.readthedocs.io/en/latest/tutorial_environments.html):
   ```bash
   spack env create mfpic-dev
   spack env activate mfpic-dev
   ```

1. Download source code:
   ```bash
   git clone git@github.com:sandialabs/mfpic.git /desired/path/to/source/
   ```

1. Set up a Spack repo to tell Spack how to build the code.
   This step is only required because the code is not in Spack's built-in package repo,
   and can be removed if/when it's added.
   ```bash
   spack repo add /path/to/source/spack_repo/mfpic
   ```
   Note that this only adds the repo for the currently activated environment.
   If you like to have multiple builds, you can set up a separate environment for each one so they don't conflict.

1. Set up the Spack environment:
   ```bash
   spack add mfpic
   spack develop --no-clone --path /absolute/path/to/source --build-directory /absolute/path/to/build/directory mfpic
   spack concretize -f --test=root
   spack install --test=root
   spack build-env mfpic -- bash
   cd /path/to/build
   ninja && ninja install
   ```

1. To return to the build environment on a fresh shell:
   ```bash
   spack env activate mfpic-dev
   spack build-env mfpic -- bash
   cd /path/to/build
   ```
   You can simply
   ```bash
   exit
   despacktivate
   ```
   to exit the build environment.

1. Once `mfpic` is built, the tests can be run with `ctest` from inside the build directory.
   End-to-end tests tend to be considerably longer than unit tests;
   they can be excluded with
   ```bash
   ctest -E end_to_end_tests
   ```

   Many end-to-end tests involve running the code and analyzing
   the output in two separate steps.
   Since the "run" step is a dependency of the "analyze" step,
   ctest will automatically add the "run" step whenever the "analyze" step is requested.
   While developing tests, however, it is sometimes useful to iterate on the analysis code
   without rerunning `mfpic` every time.
   Typically this can be done with
   ```bash
   ctest -R <test_name>:analyze -FS <test_name>
   ```

## Debug build

A debug build is supported. To make a debug build, add the `build_type=Debug` variant to `mfpic`. That is, when setting up the Spack environment (step 7), use `mfpic build_type=Debug` in place of `mfpic`.

When this variant is specified, debug builds of `mfem` and `hdf5` are used and `gdb` is available in the build environment.
