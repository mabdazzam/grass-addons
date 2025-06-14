#!/usr/bin/env python
##############################################################################
# MODULE:    r.scsrunoff
#
# AUTHOR(S): Abdullah Azzam <mabdazzam@outlook.com>
#
# PURPOSE:   Generates a runoff depth raster, total runoff volume, and peak
#            discharge raster for a rainfall event using the SCS Curve Number method
#
# COPYRIGHT: (C) 2025 by Abdullah Azzam and the GRASS Development Team
#
#            This program is free software under the GNU General Public
#            License (>=v2). Read the file COPYING that comes with GRASS
#            for details.
##############################################################################

"""generates a runoff depth raster, total runoff volume, and peak discharge raster for a rainfall event using the scs curve number method"""

# %module
# % description: Generates a runoff depth raster, total runoff volume, and peak discharge raster using the SCS Curve Number method
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
# % description: Curve Number raster (minimum 30)
# %end
# %option G_OPT_R_OUTPUT
# % key: runoff
# % description: Runoff depth raster (inches)
# %end
# %option
# % key: volume
# % type: string
# % description: Name for total runoff volume output (ft³ or m³)
# % required: no
# %end
# %option G_OPT_R_OUTPUT
# % key: peak_raster
# % type: string
# % description: Peak discharge raster (cfs or cms)
# % required: no
# %end
# %option G_OPT_R_INPUT
# % key: tc_raster
# % type: string
# % description: Time of concentration raster (hours)
# % required: no
# %end
# %option G_OPT_R_INPUT
# % key: fdr
# % type: string
# % description: Flow direction raster
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

class runoff_calculator:
    """handles scs runoff depth computation"""
    @staticmethod
    def compute_runoff_depth(rainfall_raster, cn_raster, output_raster, overwrite=False):
        """compute runoff depth using q = (p - 0.2s)^2 / (p + 0.8s), s = 1000/cn - 10"""
        rainfall = garray.array(rainfall_raster)
        cn = garray.array(cn_raster)
        valid_mask = np.isfinite(cn) & (cn > 0)
        s = np.zeros_like(cn)
        s[valid_mask] = (1000.0 / cn[valid_mask]) - 10.0
        p_minus_02s = rainfall - 0.2 * s
        denominator = rainfall + 0.8 * s
        if np.any(~valid_mask | (denominator == 0)):
            gs.warning(_("invalid curve number values (zero or null) detected. minimum cn is 30. check curve_number raster. affected cells set to null"))
        runoff = np.where(
            (p_minus_02s <= 0) | (denominator == 0) | ~np.isfinite(rainfall) | ~valid_mask,
            np.nan,
            (p_minus_02s ** 2) / denominator
        )
        output_array = garray.array()
        output_array[...] = runoff
        output_array.write(output_raster, overwrite=overwrite)

class volume_calculator:
    """computes total runoff volume"""
    @staticmethod
    def compute_runoff_volume(runoff_raster):
        """compute volume as sum of runoff depth over valid cells times cell area"""
        region = gs.region()
        cell_area = region["nsres"] * region["ewres"]
        runoff = garray.array(runoff_raster)
        valid_mask = np.isfinite(runoff)
        volume = np.sum(runoff[valid_mask] * cell_area)
        return volume if np.isfinite(volume) else 0.0

class peak_discharge_calculator:
    """computes peak discharge raster"""
    @staticmethod
    def compute_peak_discharge(runoff_raster, tc_raster, fdr_raster, duration, output_raster, overwrite=False):
        """compute peak discharge using qp = (484 * q_weighted * a) / (0.5d + 0.6tc)"""
        flow_acc = gs.append_node_pid("flow_acc")
        q_weighted = gs.append_node_pid("q_weighted")
        try:
            # compute flow accumulation (cell count)
            gs.run_command("r.flowaccumulation", input=fdr_raster, output=flow_acc, type="FCELL", overwrite=True, quiet=True)
            # compute weighted flow accumulation (q_weighted)
            gs.run_command("r.flowaccumulation", input=fdr_raster, weight=runoff_raster, output=q_weighted, type="FCELL", overwrite=True, quiet=True)
            flow_acc_array = garray.array(flow_acc)
            q_weighted_array = garray.array(q_weighted)
            tc = garray.array(tc_raster)
            region = gs.region()
            cell_area = region["nsres"] * region["ewres"] / 27878400.0
            area = flow_acc_array * cell_area
            tp = 0.5 * duration + 0.6 * tc
            valid_mask = (tp > 0) & np.isfinite(area) & np.isfinite(q_weighted_array) & np.isfinite(tp)
            qp = np.full_like(q_weighted_array, np.nan)
            qp[valid_mask] = (484.0 * q_weighted_array[valid_mask] * area[valid_mask]) / tp[valid_mask]
            output_array = garray.array()
            output_array[...] = qp
            output_array.write(output_raster, overwrite=overwrite)
        finally:
            gs.run_command("g.remove", type="raster", name=[flow_acc, q_weighted], flags="f", quiet=True)

class unit_converter:
    """handles unit conversions"""
    @staticmethod
    def convert_rainfall(rainfall_raster, output_raster):
        """convert rainfall from mm to inches"""
        gs.mapcalc(f"{output_raster} = {rainfall_raster} * {1.0 / 25.4}", overwrite=True)
        return output_raster

    @staticmethod
    def convert_volume(volume):
        """convert volume from ft³ to m³"""
        return volume * 0.0283168

    @staticmethod
    def convert_discharge_raster(input_raster, output_raster):
        """convert discharge from cfs to cms"""
        gs.mapcalc(f"{output_raster} = {input_raster} * 0.0283168", overwrite=True)

class output_writer:
    """writes output to console"""
    @staticmethod
    def write_output(output_type, value, units):
        """write volume output"""
        unit_label = "ft³" if units == "us" else "m³"
        gs.message(_("{type}: {value:.2f} {unit}").format(
            type=output_type.capitalize(), value=value, unit=unit_label))

class scs_runoff:
    """orchestrates scs runoff calculations"""
    def __init__(self, options, flags):
        self.options = options
        self.flags = flags
        self.rainfall = options["rainfall"]
        self.curve_number = options["curve_number"]
        self.runoff = options["runoff"]
        self.volume = options["volume"]
        self.peak_raster = options["peak_raster"]
        self.tc_raster = options["tc_raster"]
        self.fdr = options["fdr"]
        self.duration = float(options["duration"]) if options["duration"] else None
        self.units = options["units"]
        self.overwrite = flags["o"]
        self.temp_rasters = []

    def clean_temp_rasters(self):
        """remove temporary rasters"""
        if self.temp_rasters:
            gs.run_command("g.remove", type="raster", name=self.temp_rasters, flags="f", quiet=True)

    def validate_inputs(self):
        """validate peak discharge inputs"""
        if self.peak_raster:
            if self.tc_raster is None or self.duration is None or self.fdr is None:
                gs.fatal(_("'tc_raster', 'duration', and 'fdr' required for 'peak_raster'"))
            if self.duration <= 0:
                gs.fatal(_("duration must be positive"))

    def process_rainfall(self):
        """convert rainfall units if needed"""
        if self.units == "metric":
            temp_raster = gs.append_node_pid("rainfall_us")
            self.temp_rasters.append(temp_raster)
            atexit.register(self.clean_temp_rasters)
            gs.verbose(_("converting rainfall from mm to inches"))
            return unit_converter.convert_rainfall(self.rainfall, temp_raster)
        return self.rainfall

    def compute(self):
        """execute calculations"""
        self.validate_inputs()
        gs.use_temp_region()

        # compute runoff
        gs.verbose(_("computing runoff depth"))
        runoff_calculator.compute_runoff_depth(
            self.process_rainfall(), self.curve_number, self.runoff, self.overwrite
        )

        # compute volume
        if self.volume:
            gs.verbose(_("computing volume"))
            vol = volume_calculator.compute_runoff_volume(self.runoff)
            if self.units == "metric":
                vol = unit_converter.convert_volume(vol)
            output_writer.write_output("volume", vol, self.units)

        # compute peak discharge
        if self.peak_raster:
            gs.verbose(_("computing peak discharge"))
            peak_output = self.peak_raster
            if self.units == "metric":
                temp_peak = gs.append_node_pid("peak_cfs")
                self.temp_rasters.append(temp_peak)
                peak_output = temp_peak
            peak_discharge_calculator.compute_peak_discharge(
                self.runoff, self.tc_raster, self.fdr, self.duration, peak_output, self.overwrite
            )
            if self.units == "metric":
                gs.verbose(_("converting peak discharge to cms"))
                unit_converter.convert_discharge_raster(temp_peak, self.peak_raster)

        # save history
        gs.raster_history(self.runoff, overwrite=True)
        if self.peak_raster:
            gs.raster_history(self.peak_raster, overwrite=True)

def main():
    """main function"""
    options, flags = gs.parser()
    scs_runoff(options, flags).compute()

if __name__ == "__main__":
    sys.exit(main())
