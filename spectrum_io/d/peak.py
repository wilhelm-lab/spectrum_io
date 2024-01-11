import math

# import pyximport; pyximport.install()
# from mgf_filter.cython.mat import ceil


class Peak:
    def __init__(self, mz, intensity, delta_function, meta=[]):
        """
        mz = m/z value, intensity, left border position, delta_function(needs to react to mz)
        """
        self.mz = mz
        self.intensity = intensity
        self.left = 0.0
        self.delta_function = delta_function
        self.counts = 1
        self.meta = meta
        self.update()

    def __str__(self):
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
        """
        updates delta and calculate left border
        peaks dont have right border
        """
        # left is needed for key()
        # function will be overwritten
        self.delta = self.delta_function(self.mz)
        self.left = self.mz - self.delta
        self.ceiled_key = math.ceil(self.left)

    def key(self):
        """
        key = ceil of self border
        meaning 0.1 -> 1.0
        """
        return self.ceiled_key
