import grass.script as gs
import numpy as np

def compute_runoff_depth(rainfall_raster, cn_raster, output_raster, overwrite=False):
    """
    Compute SCS runoff depth raster using Q = (P - 0.2S)^2 / (P + 0.8S)
    where S = 1000/CN - 10.
    """
    # read input rasters
    rainfall = gs.read_raster(rainfall_raster, returnnumpy=True)
    cn = gs.read_raster(cn_raster, returnnumpy=True)

    # compute S = 1000/CN - 10
    s = np.where(cn > 0, (1000.0 / cn) - 10.0, 0.0)

    # compute runoff depth Q
    p_minus_02s = rainfall - 0.2 * s
    runoff = np.where(
        p_minus_02s <= 0,
        0.0,
        (p_minus_02s ** 2) / (rainfall + 0.8 * s)
    )

    # write output raster
    gs.write_raster(output_raster, runoff, overwrite=overwrite)
