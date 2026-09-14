"""
Suffix for development builds: '+shortsha[.dirty]', or '' if HEAD is
an exact tag match. Release builds resolve @VCS_TAG@ in update_version.py
instead, so this script's output is irrelevant for those.
"""

import subprocess


def get_suffix():
    # no vcs tag if git isn't available
    try:
        subprocess.run(["git", "status"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        return ""

    try:
        subprocess.check_output(
            ["git", "describe", "--exact-match", "--tags", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return ""
    except subprocess.CalledProcessError:
        sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .decode()
            .strip()
        )
        dirty = (
            subprocess.run(
                ["git", "diff", "--quiet", "HEAD"], capture_output=True
            ).returncode
            != 0
        )
        return f"+{sha}.dirty" if dirty else f"+{sha}"


if __name__ == "__main__":
    print(get_suffix(), end="")
