#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Example launcher for SWIFT.py with PyVista visualization enabled.
#
# Default behavior:
# - uses the built-in synthetic circular mode
# - opens a static PyVista 3D window
# - overlays reference-species velocity arrows in gold
# - writes a loopable GIF into the current directory
#
# To use a real PLUME mode file instead, uncomment SWIFT_MODE_FILE and
# update the path plus SWIFT_ROW_INDEX as needed.

#export SWIFT_MODE_FILE="/home/kgklein/Codes/PLUME_git/data/example/map_par_kpar_1_10000.mode1"
export SWIFT_MODE_FILE="/home/kgklein/Codes/PLUME/data/example/map_par_kpar_1_10000.mode1"
export SWIFT_ROW_INDEX=64

export SWIFT_N_WAVELENGTHS_X=1.0
export SWIFT_N_WAVELENGTHS_Z=1.0
export SWIFT_N_LINES=6
export SWIFT_N_POINTS_PER_LINE=32
export SWIFT_Y0=0.0

export SWIFT_CREATE_STATIC_PLOT=1
export SWIFT_CREATE_GIF=1
export SWIFT_SHOW_VELOCITY=1
export SWIFT_TITLE_KEYWORD="ICW"
export SWIFT_OUTPUT_GIF="$SCRIPT_DIR/deltaB_mode.gif"
export SWIFT_ARROW_STRIDE=8
export SWIFT_ARROW_SCALE=0.12
export SWIFT_FRAMES_PER_PERIOD=16

# Set SWIFT_SHOW_VELOCITY=0 if you want to hide the gold δU_ref arrows.
# Set SWIFT_TITLE_KEYWORD="" if you do not want a mode label like ICW in the title/filename.

python3 "$SCRIPT_DIR/SWIFT.py"
