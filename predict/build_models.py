#! /usr/bin/env python3

import sys
import os

types = ["structure", "sequence"]
for arg in sys.argv:
    if arg.startswith("--types="):
        types = [t for t in arg[8:].split(",") if t in types]
        continue