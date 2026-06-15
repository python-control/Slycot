"""cibuildwheel setup script

Usage
  python cibw_before_build.py <scipy-openblas32-version> <project-root>

Installs scipy_openblas32 version given by first command-line argument.

Copies scipy_openblas32 to build-libs in project root; project root
given by second command-line argument.
"""

import os
import shutil
import subprocess
import sys


def install_openblas32(scipy_openblas32_version):
    print(f"{__file__}: Installing scipy_openblas32=={scipy_openblas32_version}")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--progress-bar",
            "off",
            f"scipy_openblas32=={scipy_openblas32_version}",
        ],
        check=True,
    )


def copy_libraries(project_root):
    import scipy_openblas32

    src = scipy_openblas32.get_lib_dir()
    dst = os.path.join(project_root, "build-libs")
    if os.path.exists(dst):
        print(f'Deleting existing {dst}')
        shutil.rmtree(dst)

    print(f"{__file__}: Copying {src} to {dst}")
    shutil.copytree(src, dst)


def main():
    version = sys.argv[1]
    project_root = sys.argv[2]
    install_openblas32(version)
    copy_libraries(project_root)


if __name__ == "__main__":
    main()
