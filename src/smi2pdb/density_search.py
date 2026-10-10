"""Increasing-density schedules; each point is a fresh packing, not compression."""
from decimal import Decimal


def density_fractions(start, increment, maximum=1.0):
    """Include the exact upper endpoint even when the increment does not divide it."""
    start, increment, maximum = (Decimal(str(v)) for v in (start, increment, maximum))
    if not all(v.is_finite() for v in (start, increment, maximum)):
        raise ValueError("Density search values must be finite.")
    if not 0 < start <= maximum <= 1 or not 0 < increment <= 1:
        raise ValueError("Require 0 < start <= maximum <= 1 and 0 < increment <= 1.")
    fraction = start
    while True:
        yield float(fraction)
        if fraction == maximum:
            return
        fraction = min(fraction + increment, maximum)
