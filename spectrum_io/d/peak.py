import math
from typing import Any, Callable, List, Optional

# import pyximport; pyximport.install()
# from mgf_filter.cython.mat import ceil


class Peak:
    """Container class for a single fragment peak."""

    def __init__(self, mz: float, intensity: float, delta_function: Callable, meta: Optional[List[Any]] = None):
        """
        Contructor for a single fragment peak.

        A fragment peak is contructed from an mz and intensity value, as well as a delta function providing a
        function to calculate the mass tolerance window and an additional metadata list # TODO provide details
        :param mz: mass to charge ratio of peak
        :param intensity: relative intensity of the peak
        :param delta_function: callable to calculate mass tolerance window
        :param meta: metadata list
        """
        self.mz = mz
        self.intensity = intensity
        self.left = 0.0
        self.delta_function = delta_function
        self.counts = 1
        self.meta = meta if meta is not None else []
        self.update()

    def __str__(self):
        """Overrides the default str method to provide information about mz, intensity, left border, and delta."""
        return (
            "mz: "
            + str(self.mz)
            + "\n"
            + "intensity: "
            + str(self.intensity)
            + "\n"
            + "left: "
            + str(self.left)
            + "\n"
            + "delta: "
            + str(self.delta)
        )

    def update(self):
        """Updates delta and calculate left border."""
        # left is needed for key()
        # function will be overwritten
        self.delta = self.delta_function(self.mz)
        self.left = self.mz - self.delta
        self.ceiled_key = math.ceil(self.left)

    def key(self):
        """Provides access to the ceiled key, which is ceil rounded left border."""
        return self.ceiled_key
