import grass.script as gs

def write_output(output_type, value, units):
    """
    Write output to GRASS console.
    """
    unit_label = "ft³" if units == "us" else "m³" if output_type == "volume" else "cfs" if units == "us" else "cms"
    gs.message(_("{type}: {value:.2f} {unit}").format(
        type=output_type.capitalize(), value=value, unit=unit_label))

