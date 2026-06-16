"""cibuildwheel setup script

Usage
  python cibw_before_build.py <scipy-openblas32-version> <project-root>

Installs scipy_openblas32 version given by first command-line argument.

Copies scipy_openblas32 to build-libs in project root; project root
given by second command-line argument.
"""

import os
import pathlib
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


def copy_slicot_license(project_root):
    inname = os.path.join(project_root, "slycot", "src", "SLICOT-Reference", "LICENSE")
    outname = os.path.join(project_root, "LICENSE-SLICOT.txt")
    print(f'Copying license file {inname} to {outname}')
    shutil.copyfile(inname, outname)


def copy_scipy_openblas32_licenses(project_root):
    from importlib.metadata import metadata, files

    outname0 = os.path.join(project_root, "LICENSE-scipy-openblas32.txt")
    with open(outname0, "wt") as out:
        print(f'Creating license file {outname0} from scipy-openblas32 license data')
        out.write(metadata('scipy_openblas32')['License'])

    license_files = [v
                     for k, v in metadata('scipy_openblas32').items()
                     if k=='License-File']
    if len(set(pathlib.Path(p).name for p in license_files)) != len(license_files):
        # distinguish by full path?
        raise ValueError('license files with duplicate names not handled')

    ilicense = 0
    for license_file in license_files:
        matches = [p for p in files('scipy_openblas32') if license_file in str(p)]
        if len(matches) != 1:
            raise ValueError(f'Expect 1 match for {license_file}, found {len(matches)}')
        inname = matches[0].locate()
        with open(inname) as infile:
            data = infile.read()
            if data == metadata('scipy_openblas32')['License']:
                print(f'License file {inname} identical to license data, skipping')
                continue

        outname = os.path.join(project_root, f"LICENSE-scipy-openblas32-{ilicense}.txt")
        print(f'Copying license file {inname} to {outname}')
        shutil.copyfile(inname, outname)
        ilicense += 1


def main():
    version = sys.argv[1]
    project_root = sys.argv[2]
    install_openblas32(version)
    copy_libraries(project_root)
    copy_slicot_license(project_root)
    copy_scipy_openblas32_licenses(project_root)


if __name__ == "__main__":
    main()
