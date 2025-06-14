import grass.script as gs

def convert_to_us_units(rainfall_raster, output_raster):
    """
    Convert rainfall raster from mm to inches.
    """
    factor = 1.0 / 25.4  # mm to inches
    gs.mapcalc(f"{output_raster} = {rainfall_raster} * {factor}", overwrite=True)
    return output_raster

def convert_volume_to_metric(volume):
    """
    Convert volume from ft³ to m³.
    """
    return volume * 0.0283168

def convert_discharge_to_metric(discharge):
    """
    Convert discharge from cfs to cms.
    """
    return discharge * 0.0283168
