#!/usr/bin/env python3
##############################################################################
# MODULE:    r.scs.runoff
#
# AUTHOR(S): Abdullah Azzam <mabdazzam@outlook.com>
#
# PURPOSE:   Generates a runoff depth raster and optionally a peak discharge raster
#            for a rainfall event using the SCS Curve Number method, printing
#            maximum runoff volume or peak discharge values as specified
#
# COPYRIGHT: (C) 2025 by Abdullah Azzam and the GRASS Development Team
#            Licensed under the GNU General Public License (>=v2)
##############################################################################

"""
Generates a runoff depth raster and optionally a peak discharge raster for a
rainfall event using the SCS Curve Number method, printing maximum runoff
volume (cuft or cum) or peak discharge (cfs or cms) values as specified.
"""

# %module
# % description: Generates a runoff depth raster and optionally a peak discharge raster using the SCS Curve Number method
# % keyword: raster
# % keyword: hydrology
# % keyword: runoff
# %end
# %option G_OPT_R_INPUT
# % key: rainfall
# % description: Name of input rainfall depth raster (in or mm)
# %end
# %option G_OPT_R_INPUT
# % key: curve_number
# % description: Name of input curve number raster
# %end
# %option G_OPT_R_OUTPUT
# % key: runoff_depth
# % description: Name of output runoff depth raster (in or mm)
# %end
# %option
# % key: units
# % type: string
# % description: Units of analysis (us or metric)
# % options: us,metric
# % answer: us
# %end
# %option G_OPT_R_OUTPUT
# % key: peak_discharge
# % type: string
# % description: Name of output peak discharge raster (cfs or cms)
# % required: no
# %end
# %option G_OPT_R_INPUT
# % key: time_of_concentration
# % type: string
# % description: Name of input time-of-concentration raster (hrs)
# % required: no
# %end
# %option G_OPT_R_INPUT
# % key: direction
# % type: string
# % description: Name of input flow direction raster
# % required: no
# %end
# %option G_OPT_R_INPUT
# % key: accumulation
# % type: string
# % description: Name of input flow accumulation raster (cell count)
# % required: no
# %end
# %option
# % key: duration
# % type: double
# % description: Rainfall duration (hrs)
# % required: no
# %end
# %flag
# % key: o
# % description: Allow output rasters to be overwritten
# %end

import sys
import grass.script as gs

def validate_direction(dir_raster, overwrite):
    """
    Validate that the flow direction raster contains only valid D8 values (1-8).
    Returns a cleaned raster if invalid values are found.
    """
    stats = gs.parse_command("r.univar", map=dir_raster, flags="ge")
    min_val = float(stats.get("min", 0))
    max_val = float(stats.get("max", 8))
    if min_val < 1 or max_val > 8:
        gs.warning(_("Invalid flow direction values detected in {}. Values must be 1-8. Cleaning raster...").format(dir_raster))
        clean_raster = gs.append_node_pid("dir_clean")
        gs.mapcalc(
            exp=f"{clean_raster} = if({dir_raster} >= 1 && {dir_raster} <= 8, {dir_raster}, null())",
            overwrite=overwrite
        )
        return clean_raster
    return dir_raster

def compute_q_scs(rainfall_raster, cn_raster, output_raster, units, overwrite):
    """
    Compute q_scs = (p - 0.2 s)^2 / (p + 0.8 s), where s = ks * (1000 / cn - 10).
    ks = 25 for metric (mm), 1 for us (in). Output in mm (metric) or in (us).
    """
    s_raster = gs.append_node_pid("s")
    temp1 = gs.append_node_pid("p_minus")
    temp2 = gs.append_node_pid("denom")
    temp_list = [s_raster, temp1, temp2]
    try:
        ks = 25.0 if units == "metric" else 1.0
        gs.mapcalc(
            exp=(
                f"{s_raster} = if({cn_raster} >= 30,"
                f" {ks} * (1000.0/{cn_raster} - 10.0), null())"
            ),
            overwrite=overwrite
        )
        stats = gs.parse_command("r.univar", map=s_raster, flags="ge")
        if float(stats.get("null_cells", 0)) > 0 or float(stats.get("min", 0)) < 0:
            gs.warning(
                _("Invalid curve number values (<= 0 or null) detected. "
                  "Minimum CN for SCS Method (and generally in hydrology) is 30. "
                  "Affected cells set to null. Check input curve_number raster.")
            )
        gs.mapcalc(
            exp=f"{temp1} = {rainfall_raster} - 0.2 * {s_raster}",
            overwrite=overwrite
        )
        gs.mapcalc(
            exp=f"{temp2} = {rainfall_raster} + 0.8 * {s_raster}",
            overwrite=overwrite
        )
        gs.mapcalc(
            exp=(
                f"{output_raster} = if({temp1} > 0 && {temp2} != 0,"
                f" ({temp1}^2)/{temp2}, null())"
            ),
            overwrite=overwrite
        )
    finally:
        gs.run_command("g.remove", type="raster", name=temp_list, flags="f", quiet=True)

def compute_volume(q_scs_raster, units):
    """
    Compute total runoff volume = sum(q_scs * cell_area) in cum (metric) or cuft (us).
    q_scs in mm (metric) or in (us), cell_area in sqm (metric) or sqft (us).
    """
    region = gs.region()
    nsres = float(region["nsres"])
    ewres = float(region["ewres"])
    cell_area = nsres * ewres  # sqm
    if units == "us":
        cell_area *= 10.7639  # sqm to sqft
    gs.message(f"cell_area: {cell_area} {'sqm' if units == 'metric' else 'sqft'}")

    stats = gs.parse_command("r.univar", map=q_scs_raster, flags="ge")
    depth_sum = float(stats.get("sum", 0))
    volume = depth_sum * cell_area
    if units == "metric":
        volume /= 1000.0  # mm*sqm to cum
    else:
        volume /= 12.0    # in*sqft to cuft
    return volume if volume != float('inf') else 0.0

def compute_peak_discharge(q_scs_raster, tc_raster, dir_raster, fac_raster,
                           duration, output_raster, units, overwrite):
    """
    Compute peak discharge qp = (kp * q_weighted * area) / (0.5 * duration + 0.6 * tc).
    kp = 0.208 for metric (cms), 484 for us (cfs). Area in sqkm (metric) or sqmi (us).
    q_weighted in mm (metric) or in (us), tc and duration in hrs.
    """
    weighted_raster = gs.append_node_pid("qw")
    area_raster = gs.append_node_pid("area")
    tp_raster = gs.append_node_pid("tp")
    numerator = gs.append_node_pid("numerator")
    temp_list = [weighted_raster, area_raster, tp_raster, numerator]
    try:
        dir_raster = validate_direction(dir_raster, overwrite)

        gs.run_command(
            "r.flowaccumulation",
            input=dir_raster,
            weight=q_scs_raster,
            output=weighted_raster,
            type="FCELL",
            overwrite=overwrite,
            quiet=True,
        )
        qw_stats = gs.parse_command("r.univar", map=weighted_raster, flags="ge")
        gs.message(f"weighted_raster stats: min={qw_stats.get('min', 'N/A')} {'mm' if units == 'metric' else 'in'}, "
                   f"max={qw_stats.get('max', 'N/A')} {'mm' if units == 'metric' else 'in'}, "
                   f"null_cells={qw_stats.get('null_cells', 'N/A')}")

        region = gs.region()
        nsres = float(region["nsres"])
        ewres = float(region["ewres"])
        cell_area = nsres * ewres  # sqm
        area_conversion = 0.000001 if units == "metric" else 1.0 / (5280.0 * 5280.0)  # sqm to sqkm or sqmi
        gs.message(f"cell_area: {cell_area} sqm, area_conversion: {area_conversion} {'sqkm/sqm' if units == 'metric' else 'sqmi/sqm'}")

        gs.mapcalc(
            exp=f"{area_raster} = {fac_raster} * {cell_area} * {area_conversion}",
            overwrite=overwrite
        )
        area_stats = gs.parse_command("r.univar", map=area_raster, flags="ge")
        gs.message(f"area_raster stats: min={area_stats.get('min', 'N/A')} {'sqkm' if units == 'metric' else 'sqmi'}, "
                   f"max={area_stats.get('max', 'N/A')} {'sqkm' if units == 'metric' else 'sqmi'}, "
                   f"null_cells={area_stats.get('null_cells', 'N/A')}")

        # Fill null cells in tc_raster with 0.1 hrs
        tc_clean = gs.append_node_pid("tc_clean")
        gs.mapcalc(
            exp=f"{tc_clean} = if(isnull({tc_raster}), 0.1, {tc_raster})",
            overwrite=overwrite
        )
        temp_list.append(tc_clean)

        # Apply tc filter
        tc_filtered = gs.append_node_pid("tc_filtered")
        gs.mapcalc(
            exp=f"{tc_filtered} = if({tc_clean} >= 0.1, {tc_clean}, 0.1)",
            overwrite=overwrite
        )
        temp_list.append(tc_filtered)

        gs.mapcalc(
            exp=f"{tp_raster} = 0.5 * {duration} + 0.6 * {tc_filtered}",
            overwrite=overwrite
        )
        tc_stats = gs.parse_command("r.univar", map=tc_filtered, flags="ge")
        tp_stats = gs.parse_command("r.univar", map=tp_raster, flags="ge")
        gs.message(f"tc_filtered stats: min={tc_stats.get('min', 'N/A')} hrs, "
                   f"max={tc_stats.get('max', 'N/A')} hrs, "
                   f"null_cells={tc_stats.get('null_cells', 'N/A')}")
        gs.message(f"tp_raster stats: min={tp_stats.get('min', 'N/A')} hrs, "
                   f"max={tp_stats.get('max', 'N/A')} hrs, "
                   f"null_cells={tp_stats.get('null_cells', 'N/A')}")

        # Tighter numerator cap to prevent overflow
        gs.mapcalc(
            exp=f"{numerator} = min({weighted_raster} * {area_raster}, 100000)",
            overwrite=overwrite
        )
        num_stats = gs.parse_command("r.univar", map=numerator, flags="ge")
        gs.message(f"numerator stats: min={num_stats.get('min', 'N/A')} {'mm*sqkm' if units == 'metric' else 'in*sqmi'}, "
                   f"max={num_stats.get('max', 'N/A')} {'mm*sqkm' if units == 'metric' else 'in*sqmi'}")

        kp = 0.208 if units == "metric" else 484.0
        gs.mapcalc(
            exp=(
                f"{output_raster} = if({tp_raster} > 0.1 && {weighted_raster} != null() "
                f"&& {area_raster} != null(), ({kp} * {numerator}) / {tp_raster}, null())"
            ),
            overwrite=overwrite
        )
    finally:
        gs.run_command("g.remove", type="raster", name=temp_list, flags="f", quiet=True)

def log_max_values(q_scs_raster, peak_raster=None, units="us"):
    """
    Log maximum runoff volume (cum or cuft) or peak discharge (cms or cfs).
    """
    if peak_raster:
        stats = gs.parse_command("r.univar", map=peak_raster, flags="ge")
        max_val = float(stats.get("max", 0))
        unit_label = "cfs" if units == "us" else "cms"
        gs.message(f"maximum peak discharge: {max_val:.2f} {unit_label}")
    volume = compute_volume(q_scs_raster, units)
    unit_label = "cuft" if units == "us" else "cum"
    gs.message(f"maximum runoff volume: {volume:.2f} {unit_label}")

class scs_runoff:
    def __init__(self, options, flags):
        self.opts = options
        self.flg = flags
        self.rainfall = options['rainfall']
        self.curve_number = options['curve_number']
        self.runoff_depth = options['runoff_depth']
        self.units = options['units']
        self.peak_discharge = options.get('peak_discharge')
        self.time_of_concentration = options.get('time_of_concentration')
        self.direction = options.get('direction')
        self.accumulation = options.get('accumulation')
        self.duration = float(options['duration']) if options['duration'] else None
        self.overwrite = flags['o']

    def validate_inputs(self):
        if self.peak_discharge:
            req = [self.time_of_concentration, self.duration,
                   self.direction, self.accumulation]
            if any(x is None for x in req):
                gs.fatal(
                    "time_of_concentration, duration, direction, and accumulation "
                    "are required for peak_discharge"
                )
            if self.duration <= 0:
                gs.fatal("duration must be positive")

    def compute(self):
        self.validate_inputs()
        gs.use_temp_region()

        gs.verbose("computing runoff depth (q_scs)")
        compute_q_scs(
            self.rainfall, self.curve_number,
            self.runoff_depth, self.units, self.overwrite
        )

        gs.verbose("logging maximum runoff volume")
        log_max_values(self.runoff_depth, units=self.units)

        if self.peak_discharge:
            gs.verbose("computing peak discharge (qp)")
            compute_peak_discharge(
                self.runoff_depth, self.time_of_concentration,
                self.direction, self.accumulation,
                self.duration, self.peak_discharge,
                self.units, self.overwrite
            )
            gs.verbose("logging maximum peak discharge")
            log_max_values(self.runoff_depth, self.peak_discharge, self.units)

        gs.raster_history(self.runoff_depth, overwrite=self.overwrite)
        if self.peak_discharge:
            gs.raster_history(self.peak_discharge, overwrite=self.overwrite)

def main():
    options, flags = gs.parser()
    scs_runoff(options, flags).compute()

if __name__ == '__main__':
    main()
