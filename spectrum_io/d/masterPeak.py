from math import ceil
from typing import TypeVar

from .peak import Peak

# import pyximport; pyximport.install()
# from mgf_filter.cython.mat import ceil


MasterPeakT = TypeVar("MasterPeakT", bound="MasterPeak")


class MasterPeak(Peak):
    """Container class for a aggregated, summed up fragment peak."""

    def __init__(self, peak: Peak):
        """
        Contructor for a master peak.

        # TODO provide details

        :param peak: a peak object # TODO provide details
        """
        # right and counts are calculated based on update func
        self.right = 0
        self.counts = 1
        self.rel_intensity_ratio = 0.0
        # ratio is set to an actual value first time by comparing a spectrum to another spectrum
        self.counts_ratio = 0.0
        self.mz_origin = peak.mz
        super().__init__(peak.mz, peak.intensity, peak.delta_function, meta=peak.meta)
        # update does not have to be called, because constructor
        # of peak calls update of MasterPeak
        # self.update()

    def __eq__(self, other: MasterPeakT) -> bool:  # type: ignore
        """
        Reports true if both master peaks have the same member variables!

        :param other: the other master peak object to compare this object with
        :return: whether or not the two objects are equal
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

    def __ne__(self, other: MasterPeakT) -> bool:  # type: ignore
        """
        Reports true if any member variable differs between the two master peaks!

        :param other: the other master peak object to compare this object with
        :return: whether or not the two objects are equal
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
            return True  # TODO check if this is really true, dicts equal means the objects should be equal

    def __str__(self) -> str:
        """
        Overrides the default str method to provide information about member variables.

        :return: string representation of this object
        """
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
        """Calculates delta, left and right."""
        self.delta = self.delta_function(self.mz)
        self.left = self.mz - self.delta
        self.right = self.mz + self.delta
        self.ceiled_key = ceil(self.left)

    def is_inside(self, peak: Peak) -> bool:
        """
        Reports true if tested peak is inside window.

        :param peak: a peak object to check

        :return: whether the peak's mz is within left and right border
        """
        return (self.left < peak.mz) and (self.right > peak.mz)

    def is_inside_mz(self, mz: float) -> bool:
        """
        Reports true if tested mz is inside window.

        :param mz: the mz to check for

        :return: whether mz is within left and right border
        """
        return (self.left < mz) and (self.right > mz)

    def add(self, peak: Peak):
        """
        Default_counts must be overwritten if more than one peak is added (multimerge).

        :param peak: a peak object to be added to this master peak.
        """
        self.mz = (self.mz * self.intensity + peak.mz * peak.intensity) / (peak.intensity + self.intensity)
        self.intensity = self.intensity + peak.intensity
        self.counts += peak.counts
        self.update()

    def smaller(self, peak: Peak) -> bool:
        """
        Report true if MasterPeak window is below tested peak  # TODO make this a magic function.

        :param peak: the other peak to compare to
        :return: wheter or not the peak is smaller to this object
        """
        return self.right < peak.mz

    def greater(self, peak: Peak):
        """
        Report true if MasterPeak window is above tested peak  # TODO make this a magic function.

        :param peak: the other peak to compare to
        :return: wheter or not the peak is greater to this object
        """
        return self.left > peak.mz

    def recalculate_ratio(self, mp: MasterPeakT):
        """
        Recalculate the ratio  # TODO details.

        :param mp: the other master peak object
        """
        rel_ratio = (self.intensity / self.counts) / (mp.intensity / mp.counts)
        self.rel_intensity_ratio = rel_ratio
        rel_counts_ratio = self.counts / mp.counts
        self.counts_ratio = rel_counts_ratio
