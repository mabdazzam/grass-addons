import grass.script as gs
import numpy as np

def compute_peak_discharge(runoff_raster, tc, duration):
    """
    Compute peak discharge using SCS triangular unit hydrograph (in cfs).
    """
    # get watershed area
    region = gs.region()
    area_sqmi = (region["rows"] * region["cols"] * region["nsres"] * region["ewres"]) / 27878400.0  # ft² to mi²

    # compute average runoff depth
    runoff = gs.read_raster(runoff_raster, returnnumpy=True)
    avg_runoff = np.mean(runoff[runoff > 0]) if np.any(runoff > 0) else 0.0

    # compute time to peak
    tp = 0.5 * duration + 0.6 * tc

    # compute peak discharge (US units: cfs)
    qp = 0.0
    if tp > 0:
        qp = (484.0 * area_sqmi * avg_runoff) / tp
    return qp
