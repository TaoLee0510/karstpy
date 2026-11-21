"""Verify the macOS 15 platform tag and native extension in each CI wheel."""

import argparse
from email import message_from_bytes
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile


ARCHITECTURES = {"ARM64": "arm64", "X64": "x86_64"}
DEPLOYMENT_TARGET = "15.0"


def check_wheel(wheel: Path, architecture: str) -> None:
    platform_tag = f"macosx_{DEPLOYMENT_TARGET.replace('.', '_')}_{architecture}"
    _, python_tag, abi_tag, filename_platform = wheel.stem.rsplit("-", 3)
    if filename_platform != platform_tag:
        raise ValueError(f"expected filename platform tag {platform_tag}, found {filename_platform}")
    expected_tag = f"{python_tag}-{abi_tag}-{platform_tag}"

    with ZipFile(wheel) as archive:
        metadata, = (name for name in archive.namelist() if name.endswith(".dist-info/WHEEL"))
        tags = message_from_bytes(archive.read(metadata)).get_all("Tag", [])
        if tags != [expected_tag]:
            raise ValueError(f"expected WHEEL tag {expected_tag}, found {tags}")
        extension, = (name for name in archive.namelist()
                      if name.startswith("karstpy/_core.") and name.endswith(".so"))
        with TemporaryDirectory() as directory:
            binary = Path(directory) / Path(extension).name
            binary.write_bytes(archive.read(extension))
            arch = subprocess.check_output(["lipo", "-archs", binary], text=True).strip()
            if arch != architecture:
                raise ValueError(f"expected Mach-O architecture {architecture}, found {arch}")
            build = subprocess.check_output(["xcrun", "vtool", "-show-build", binary], text=True)

    platforms = re.findall(r"^\s*platform\s+(\S+)\s*$", build, re.MULTILINE)
    targets = re.findall(r"^\s*minos\s+(\S+)\s*$", build, re.MULTILINE)
    if platforms != ["MACOS"] or targets != [DEPLOYMENT_TARGET]:
        raise ValueError(f"expected Mach-O minimum macOS {DEPLOYMENT_TARGET}\n{build}")
    print(f"{wheel.name}: {expected_tag}, Mach-O minimum macOS {DEPLOYMENT_TARGET}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheelhouse", type=Path)
    parser.add_argument("runner_arch", choices=ARCHITECTURES)
    args = parser.parse_args()
    wheels = sorted(args.wheelhouse.glob("*.whl"))
    if not wheels:
        raise SystemExit(f"No wheels found in {args.wheelhouse}")
    for wheel in wheels:
        try:
            check_wheel(wheel, ARCHITECTURES[args.runner_arch])
        except (ValueError, subprocess.CalledProcessError) as error:
            raise SystemExit(f"{wheel.name}: {error}") from error


if __name__ == "__main__":
    main()
