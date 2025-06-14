import grass.script as gs
import numpy as np

def compute_runoff_volume(runoff_raster):
    """
    Compute total runoff volume by integrating runoff depth over watershed area (in ft³).
    """
    # get region settings for cell area
    region = gs.region()
    cell_area = region["nsres"] * region["ewres"]  # area in ft²

    # read runoff depth raster
    runoff = gs.read_raster(runoff_raster, returnnumpy=True)

    # compute volume: sum(runoff_depth * cell_area)
    volume = np.sum(runoff * cell_area)
    return volume
