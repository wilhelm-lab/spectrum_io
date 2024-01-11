from math import ceil

from .peak import Peak

# import pyximport; pyximport.install()
# from mgf_filter.cython.mat import ceil


class MasterPeak(Peak):
    def __init__(self, Peak):
        """
        mz, intensity, delta_function by peak
        ratio is set to an actual value first time by comparing a spectrum to another spectrum
        right and counts are calculated based on update func
        """
        self.right = 0
        self.counts = 1
        self.rel_intensity_ratio = 0
        self.counts_ratio = 0
        self.mz_origin = Peak.mz
        super().__init__(Peak.mz, Peak.intensity, Peak.delta_function, meta=Peak.meta)
        # update does not have to be called, because constructor
        # of peak calls update of MasterPeak
        # self.update()

    def __eq__(self, other):
        """
        reports true if both have the same member variables!
        normally you test for the same memory address
        """
        if self.__dict__ == other.__dict__:
            return True
        else:
            if (
                (self.counts == other.counts)
                & (self.mz == other.mz)
                & (self.left == other.left)
                & (self.right == other.right)
                & (self.mz_origin == other.mz_origin)
            ):
                return True
            else:
                return False

    def __ne__(self, other):
        """
        because we override rich comparison operator, we have to implement the other ones too
        """
        if self.__dict__ != other.__dict__:
            if (
                (self.counts == other.counts)
                & (self.mz == other.mz)
                & (self.left == other.left)
                & (self.right == other.right)
                & (self.mz_origin == other.mz_origin)
            ):
                return False
            else:
                return True
        else:
            return True

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
            + "right: "
            + str(self.right)
            + "\n"
            + "counts: "
            + str(self.counts)
            + "\n"
            + "rel_intensity_ratio: "
            + str(self.rel_intensity_ratio)
            + "\n"
            + "counts_ratio: "
            + str(self.counts_ratio)
            + "\n"
            + "origin: "
            + str(self.mz_origin)
        )

    def update(self):
        """
        calculates delta, left and right
        """
        self.delta = self.delta_function(self.mz)
        self.left = self.mz - self.delta
        self.right = self.mz + self.delta
        self.ceiled_key = ceil(self.left)

    def isInside(self, peak):
        """
        reports true if tested peak is inside window
        """
        return (self.left < peak.mz) and (self.right > peak.mz)

    def isInsideMz(self, mz):
        """
        reports true if tested my is inside window
        """
        return (self.left < mz) and (self.right > mz)

    def add(self, peak):
        """
        default_counts must be overwritten
        if more than one peak is added (multimerge)
        """
        self.mz = (self.mz * self.intensity + peak.mz * peak.intensity) / (peak.intensity + self.intensity)
        self.intensity = self.intensity + peak.intensity
        self.counts += peak.counts
        self.update()

    def smaller(self, peak):
        """
        reports true if MasterPeak window is below tested peak
        """
        return self.right < peak.mz

    def greater(self, peak):
        """
        reports true if MasterPeak window is above tested peak
        """
        return self.left > peak.mz

    def recalculate_ratio(self, mp):
        rel_ratio = (self.intensity / self.counts) / (mp.intensity / mp.counts)
        self.rel_intensity_ratio = rel_ratio
        rel_counts_ratio = self.counts / mp.counts
        self.counts_ratio = rel_counts_ratio
