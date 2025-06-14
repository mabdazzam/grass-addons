#!/usr/bin/env python
##############################################################################
# MODULE:    r.scsrunoff
#
# AUTHOR(S): Abdullah Azzam <mabdazzam@outlook.com>
#
# PURPOSE:   Generates a runoff depth raster and estimates total runoff volume
#            and peak discharge for a given rainfall event over a watershed
#
# COPYRIGHT: (C) 2025 by Abdullah Azzam and the GRASS Development Team
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
import grass.script.array as garray
import numpy as np

class RunoffCalculator:
    """Handles computation of SCS runoff depth raster."""
    @staticmethod
    def validate_curve_number(cn_raster):
        """Validate curve number raster to ensure values are in valid range."""
        cn = garray.array(cn_raster)
        if np.any(cn <= 0):
            gs.fatal(_("Curve number raster contains zero or negative values, which are invalid"))
        if np.any(cn > 100):
            gs.warning(_("Curve number raster contains values > 100, which may produce unexpected results"))

    @staticmethod
    def compute_runoff_depth(rainfall_raster, cn_raster, output_raster, overwrite=False):
        """
        Compute SCS runoff depth using Q = (P - 0.2S)^2 / (P + 0.8S)
        where S = 1000/CN - 10.
        """
        # Validate curve number
        RunoffCalculator.validate_curve_number(cn_raster)
        # Read rasters into NumPy arrays
        rainfall = garray.array(rainfall_raster)
        cn = garray.array(cn_raster)
        # Compute storage (S)
        s = np.where(cn > 0, (1000.0 / cn) - 10.0, 0.0)
        p_minus_02s = rainfall - 0.2 * s
        # Compute runoff depth, handling invalid cases
        denominator = rainfall + 0.8 * s
        runoff = np.where(
            (p_minus_02s <= 0) | (denominator == 0),
            0.0,
            (p_minus_02s ** 2) / denominator
        )
        # Replace NaN/infinite values with 0
        runoff = np.where(np.isfinite(runoff), runoff, 0.0)
        # Write output raster
        output_array = garray.array()
        output_array[...] = runoff
        output_array.write(output_raster, overwrite=overwrite)

class VolumeCalculator:
    """Computes total runoff volume over the watershed."""
    @staticmethod
    def compute_runoff_volume(runoff_raster):
        """Compute volume by integrating runoff depth over area (in ft³)."""
        region = gs.region()
        cell_area = region["nsres"] * region["ewres"]  # Area in ft²
        runoff = garray.array(runoff_raster)
        volume = np.sum(runoff * cell_area)
        return volume if np.isfinite(volume) else 0.0

class PeakDischargeCalculator:
    """Estimates peak discharge using SCS triangular unit hydrograph."""
    @staticmethod
    def compute_peak_discharge(runoff_raster, tc, duration):
        """Compute peak discharge (in cfs)."""
        region = gs.region()
        area_sqmi = (region["rows"] * region["cols"] * region["nsres"] * region["ewres"]) / 27878400.0  # ft² to mi²
        runoff = garray.array(runoff_raster)
        avg_runoff = np.mean(runoff[runoff > 0]) if np.any(runoff > 0) else 0.0
        tp = 0.5 * duration + 0.6 * tc
        qp = 0.0
        if tp > 0:
            qp = (484.0 * area_sqmi * avg_runoff) / tp
        return qp if np.isfinite(qp) else 0.0

class UnitConverter:
    """Handles unit conversions between US and metric systems."""
    @staticmethod
    def convert_to_us_units(rainfall_raster, output_raster):
        """Convert rainfall raster from mm to inches."""
        factor = 1.0 / 25.4  # mm to inches
        gs.mapcalc(f"{output_raster} = {rainfall_raster} * {factor}", overwrite=True)
        return output_raster

    @staticmethod
    def convert_volume_to_metric(volume):
        """Convert volume from ft³ to m³."""
        return volume * 0.0283168

    @staticmethod
    def convert_discharge_to_metric(discharge):
        """Convert discharge from cfs to cms."""
        return discharge * 0.0283168

class OutputWriter:
    """Manages writing of output results to GRASS console."""
    @staticmethod
    def write_output(output_type, value, units):
        """Write output to GRASS console."""
        unit_label = "ft³" if units == "us" else "m³" if output_type == "volume" else "cfs" if units == "us" else "cms"
        gs.message(_("{type}: {value:.2f} {unit}").format(
            type=output_type.capitalize(), value=value, unit=unit_label))

class SCSSRunoff:
    """Main class orchestrating the SCS runoff calculations."""
    def __init__(self, options, flags):
        self.options = options
        self.flags = flags
        self.rainfall = options["rainfall"]
        self.curve_number = options["curve_number"]
        self.runoff = options["runoff"]
        self.volume = options["volume"]
        self.peak = options["peak"]
        self.tc = float(options["tc"]) if options["tc"] else None
        self.duration = float(options["duration"]) if options["duration"] else None
        self.units = options["units"]
        self.overwrite = flags["o"]
        self.temp_raster = None

    def clean_temp_raster(self):
        """Remove temporary raster."""
        if self.temp_raster:
            gs.run_command("g.remove", type="raster", name=self.temp_raster, flags="f", quiet=True)

    def validate_inputs(self):
        """Validate inputs for peak discharge."""
        if self.peak and (self.tc is None or self.duration is None):
            gs.fatal(_("Both 'tc' and 'duration' are required when 'peak' is specified"))
        if self.tc is not None and self.tc <= 0:
            gs.fatal(_("Time of concentration (tc) must be positive"))
        if self.duration is not None and self.duration <= 0:
            gs.fatal(_("Rainfall duration must be positive"))

    def process_rainfall(self):
        """Handle rainfall unit conversion if needed."""
        rainfall_input = self.rainfall
        if self.units == "metric":
            self.temp_raster = gs.append_node_pid("rainfall_us")
            atexit.register(self.clean_temp_raster)
            gs.verbose(_("Converting rainfall from mm to inches"))
            rainfall_input = UnitConverter.convert_to_us_units(self.rainfall, self.temp_raster)
        return rainfall_input

    def compute(self):
        """Execute the SCS runoff calculations."""
        self.validate_inputs()
        gs.use_temp_region()

        # Compute runoff depth
        gs.verbose(_("Computing runoff depth raster"))
        RunoffCalculator.compute_runoff_depth(
            self.process_rainfall(), self.curve_number, self.runoff, overwrite=self.overwrite
        )

        # Compute total runoff volume
        if self.volume:
            gs.verbose(_("Computing total runoff volume"))
            vol = VolumeCalculator.compute_runoff_volume(self.runoff)
            if self.units == "metric":
                vol = UnitConverter.convert_volume_to_metric(vol)
            OutputWriter.write_output("volume", vol, self.units)

        # Compute peak discharge
        if self.peak:
            gs.verbose(_("Computing peak discharge"))
            discharge = PeakDischargeCalculator.compute_peak_discharge(self.runoff, self.tc, self.duration)
            if self.units == "metric":
                discharge = UnitConverter.convert_discharge_to_metric(discharge)
            OutputWriter.write_output("peak_discharge", discharge, self.units)

        # Save history
        gs.raster_history(self.runoff, overwrite=True)

def main():
    """Main function to run the SCS runoff module."""
    options, flags = gs.parser()
    scs_runoff = SCSSRunoff(options, flags)
    scs_runoff.compute()

if __name__ == "__main__":
    sys.exit(main())
