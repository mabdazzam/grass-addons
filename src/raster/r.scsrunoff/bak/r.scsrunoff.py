#!/usr/bin/env python

##############################################################################
# mODULE:    r.scsrunoff
#
# aUTHOR(S): Abdullah Azzam <mabdazzam@outlook.com>
#
# pURPOSE:   Generates a runoff depth raster and estimates total runoff volume and peak discharge for a given rainfall event over a watershed
#
# cOPYRIGHT: (C) 2025 by Abdullah Azzam and the GRASS Development Team
#
#            This program is free software under the GNU General Public
#            License (>=v2). Read the file COPYING that comes with GRASS
#            for details.
##############################################################################

"""Generates a runoff depth raster and estimates total runoff volume and peak discharge for a given rainfall event over a watershed"""

# %module
# % description: Generates a runoff depth raster and estimates total runoff volume and peak discharge using the SCS Curve Number method
# % keyword: raster
# % keyword: hydrology
# % keyword: runoff
# %end
# %option G_OPT_R_INPUT
# % key: rainfall
# % description: Rainfall depth raster (inches or mm)
# %end
# %option G_OPT_R_INPUT
# % key: curve_number
# % description: Curve Number raster
# %end
# %option G_OPT_R_OUTPUT
# % key: runoff
# % description: Runoff depth raster
# %end
# %option
# % key: volume
# % type: string
# % description: Name for total runoff volume output (ft³ or m³)
# % required: no
# %end
# %option
# % key: peak
# % type: string
# % description: Name for peak discharge output (cfs or cms)
# % required: no
# %end
# %option
# % key: tc
# % type: double
# % description: Time of concentration (hours)
# % required: no
# %end
# %option
# % key: duration
# % type: double
# % description: Rainfall duration (hours)
# % required: no
# %end
# %option
# % key: units
# % type: string
# % description: Output units
# % options: us,metric
# % answer: us
# %end
# %flag
# % key: o
# % description: Overwrite existing output rasters
# %end

import sys
import atexit
import grass.script as gs
from grass.script import parser, run_command, fatal
from runoff_calc import compute_runoff_depth
from volume_calc import compute_runoff_volume
from peak_discharge import compute_peak_discharge
from utils import convert_to_us_units, convert_volume_to_metric, convert_discharge_to_metric
from io_utils import write_output

def clean(name):
    """Remove temporary raster."""
    gs.run_command("g.remove", type="raster", name=name, flags="f", quiet=True)

def main():
    """Main function to orchestrate SCS runoff calculations."""
    options, flags = parser()
    rainfall = options["rainfall"]
    curve_number = options["curve_number"]
    runoff = options["runoff"]
    volume = options["volume"]
    peak = options["peak"]
    tc = float(options["tc"]) if options["tc"] else None
    duration = float(options["duration"]) if options["duration"] else None
    units = options["units"]
    overwrite = flags["o"]

    # validate inputs for peak discharge
    if peak and (tc is None or duration is None):
        fatal(_("Both 'tc' and 'duration' are required when 'peak' is specified"))

    # create temporary raster for unit conversion
    temp_raster = None
    if units == "metric":
        temp_raster = gs.append_node_pid("rainfall_us")
        atexit.register(clean, temp_raster)

    gs.use_temp_region()

    # convert rainfall to US units if needed
    rainfall_input = rainfall
    if units == "metric":
        gs.verbose(_("Converting rainfall from mm to inches"))
        rainfall_input = convert_to_us_units(rainfall, temp_raster)

    # compute runoff depth
    gs.verbose(_("Computing runoff depth raster"))
    compute_runoff_depth(rainfall_input, curve_number, runoff, overwrite=overwrite)

    # compute total runoff volume if requested
    if volume:
        gs.verbose(_("Computing total runoff volume"))
        vol = compute_runoff_volume(runoff)
        if units == "metric":
            vol = convert_volume_to_metric(vol)
        write_output("volume", vol, units)

    # compute peak discharge if requested
    if peak:
        gs.verbose(_("Computing peak discharge"))
        discharge = compute_peak_discharge(runoff, tc, duration)
        if units == "metric":
            discharge = convert_discharge_to_metric(discharge)
        write_output("peak_discharge", discharge, units)

    # save history
    gs.raster_history(runoff, overwrite=True)

if __name__ == "__main__":
    sys.exit(main())
