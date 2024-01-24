import csv
from math import ceil
from typing import Callable, Optional

from sortedcontainers import SortedList

from .masterPeak import MasterPeak
from .peak import Peak


def _calculate_delta_by_ppm(ppm):
    def fix_ppm(mz):
        return ppm * float(mz) / (pow(10, 6))

    return fix_ppm


def _calculate_relative_intensity(a_intensity):
    max_intensity = max(a_intensity)
    return [x / max_intensity for x in a_intensity]


class MasterSpectrum:
    """
    Container class for a Spectrum.

    Master spectrum takes peaks and converts them to master peaks. If there are no master peaks
    to be inserted spectrum is {charge of spectrum: {bins: MP}}.
    """

    def __init__(self):
        """
        Constructor for MasterSpectrum.

        :param spectrum: TODO
        :param ignore_charges: TODO
        :param delta_func: TODO
        """
        self.spectrum = {}
        self.merged = 0
        self.appended = 0
        self.multimerged = 0  # within one insertion 3 peaks are merged

    def binary(self, peak, imin, imax, charge):
        """
        Input values:
        peak
        imin is minimum search position
        imax is max search poistion
        charge defines in which masterspectrum to search
        Return values:
        first argument: position of insert, -1 if can not be added to any peak
        sc argument : add peak from left or right bin
        rd argument : must left peak also be added (merge case)
        4th argument: right peak also be added (think about 3 peaks and a merge between 2 and 3)
        """
        key = peak.key()
        imid = int(ceil((imax + imin) / 2))
        if imax < imin:
            exist_left_bin = key - 1 in self.spectrum[charge].keys()
            exist_right_bin = key + 1 in self.spectrum[charge].keys()

            if exist_left_bin:
                if self.spectrum[charge][key - 1][-1].isInside(peak):
                    return -1, -1, False, False

            if exist_right_bin:
                if self.spectrum[charge][key + 1][0].isInside(peak):
                    return -1, 1, False, False

            return -1, 0, False, False

        # search must go on
        mspeak_in_bin = self.spectrum[charge][key][imid]
        if mspeak_in_bin.greater(peak):
            return self.binary(peak, imin, imid - 1, charge)
        elif mspeak_in_bin.smaller(peak):
            return self.binary(peak, imid + 1, imax, charge)
        # search results in peak that should be added
        else:
            if imid == 0:
                exist_left_bin = key - 1 in self.spectrum[charge]
                exist_right_bin = key + 1 in self.spectrum[charge]
                if exist_left_bin:
                    if self.spectrum[charge][key - 1][-1].isInside(peak):
                        return 0, -1, False, False
                    else:
                        return 0, 0, False, False  # peak must be added to peak in pos1 but but left bin can be ignored
                elif exist_right_bin:
                    if self.spectrum[charge][key + 1][0].isInside(peak):
                        return 0, 1, False, False
                    else:
                        return 0, 0, False, False  # peak must be added to peak in pos1 but but left bin can be ignored
                else:
                    return 0, 0, False, False

            else:  # peak is somewhere between 1 and last
                # imid -1 exists alway
                is_last_entry = len(self.spectrum[charge][key]) - 1 == imid
                exists_bigger_peak = not (is_last_entry)
                if exists_bigger_peak:
                    if self.spectrum[charge][key][imid - 1].isInside(peak):
                        return imid, 0, True, False
                    if self.spectrum[charge][key][imid + 1].isInside(peak):
                        return imid, 0, False, True
                    else:
                        return imid, 0, False, False
                else:  # is last entry
                    exist_right_bin = key + 1 in self.spectrum[charge]
                    if self.spectrum[charge][key][imid - 1].isInside(peak):
                        return imid, 0, True, False
                    elif exist_right_bin:
                        if self.spectrum[charge][key + 1][0].isInside(peak):
                            return imid, 1, False, False
                        else:
                            return imid, 0, False, False
                    else:
                        return imid, 0, False, False

    def add(self, peak, charge: int = 0):
        """
        Add charges to peaks.

        Charge is by default 0, so if a spectrum should be summed without including charge information
        it is defaulted to 0.

        :param peak: TODO
        :param charge: TODO
        """
        key = peak.key()
        if charge not in self.spectrum:
            self.spectrum[charge] = {}

        if key in self.spectrum[charge]:
            idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = self.binary(
                peak, 0, len(self.spectrum[charge][key]) - 1, charge
            )
            if idx == -1:  # does not have to react to merge cases!
                if bin_to_ack == 0:
                    self.appended += 1
                    if type(peak).__name__ == "Peak":
                        self.spectrum[charge][key].add(MasterPeak(peak))
                    else:
                        self.spectrum[charge][key].add(peak)
                elif bin_to_ack == -1:
                    if len(self.spectrum[charge][key]) == 0:
                        del self.spectrum[charge][key]
                    get_master_peak_left_bin = self.spectrum[charge][key - 1][-1]
                    del self.spectrum[charge][key - 1][-1]
                    if len(self.spectrum[charge][key - 1]) == 0:
                        del self.spectrum[charge][key - 1]
                    get_master_peak_left_bin.add(peak)
                    self.merged += 1
                    if get_master_peak_left_bin.key() not in self.spectrum[charge]:
                        self.spectrum[charge][get_master_peak_left_bin.key()] = SortedList(key=lambda i: i.left)
                    self.spectrum[charge][get_master_peak_left_bin.key()].add(get_master_peak_left_bin)
                else:  # +1
                    if len(self.spectrum[charge][key]) == 0:
                        del self.spectrum[charge][key]
                    get_master_peak_right_bin = self.spectrum[charge][key + 1][0]
                    del self.spectrum[charge][key + 1][0]
                    if len(self.spectrum[charge][key + 1]) == 0:
                        del self.spectrum[charge][key + 1]
                    get_master_peak_right_bin.add(peak)
                    self.merged += 1
                    if get_master_peak_right_bin.key() not in self.spectrum[charge]:
                        self.spectrum[charge][get_master_peak_right_bin.key()] = SortedList(key=lambda i: i.left)
                    self.spectrum[charge][get_master_peak_right_bin.key()].add(get_master_peak_right_bin)

            else:  # idx != -1
                if bin_to_ack == 0:
                    if should_merge_left_peak:
                        get_left_master_peak = self.spectrum[charge][key][idx - 1]
                        get_idx_master_peak = self.spectrum[charge][key][idx]
                        del self.spectrum[charge][key][idx - 1 : idx + 1]  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del self.spectrum[charge][key]
                        get_idx_master_peak.add(get_left_master_peak)
                        get_idx_master_peak.add(peak)
                        self.multimerged += 1
                        if get_idx_master_peak.key() not in self.spectrum[charge]:
                            self.spectrum[charge][get_idx_master_peak.key()] = SortedList(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_master_peak.key()].add(get_idx_master_peak)
                    elif should_merge_right_peak:
                        get_idx_master_peak = self.spectrum[charge][key][idx]
                        get_right_master_peak = self.spectrum[charge][key][idx + 1]
                        del self.spectrum[charge][key][idx : idx + 2]  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del self.spectrum[charge][key]
                        get_idx_master_peak.add(get_right_master_peak)
                        get_idx_master_peak.add(peak)
                        self.multimerged += 1
                        if get_idx_master_peak.key() not in self.spectrum[charge]:
                            self.spectrum[charge][get_idx_master_peak.key()] = SortedList(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_master_peak.key()].add(get_idx_master_peak)
                    else:  # case idx = 0, 0, False, False
                        get_idx_master_peak = self.spectrum[charge][key][idx]
                        del self.spectrum[charge][key][idx]  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del self.spectrum[charge][key]
                        get_idx_master_peak.add(peak)
                        self.merged += 1
                        if get_idx_master_peak.key() not in self.spectrum[charge]:
                            self.spectrum[charge][get_idx_master_peak.key()] = SortedList(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_master_peak.key()].add(get_idx_master_peak)
                elif bin_to_ack == -1:  # idx != -1
                    # 0, -1, F, F
                    get_idx_master_peak = self.spectrum[charge][key][idx]
                    get_master_peak_before = self.spectrum[charge][key - 1][-1]
                    del self.spectrum[charge][key][idx]
                    if len(self.spectrum[charge][key]) == 0:
                        del self.spectrum[charge][key]
                    del self.spectrum[charge][key - 1][-1]
                    if len(self.spectrum[charge][key - 1]) == 0:
                        del self.spectrum[charge][key - 1]
                    get_idx_master_peak.add(get_master_peak_before)
                    get_idx_master_peak.add(peak)
                    self.merged += 1
                    if get_idx_master_peak.key() not in self.spectrum[charge]:
                        self.spectrum[charge][get_idx_master_peak.key()] = SortedList(key=lambda i: i.left)
                    self.spectrum[charge][get_idx_master_peak.key()].add(get_idx_master_peak)
                else:  # bin_to_ack == 1, idx !=  -1
                    get_idx_master_peak = self.spectrum[charge][key][idx]
                    get_master_peak_after = self.spectrum[charge][key + 1][0]
                    del self.spectrum[charge][key][idx]
                    if len(self.spectrum[charge][key]) == 0:
                        del self.spectrum[charge][key]
                    del self.spectrum[charge][key + 1][0]
                    if len(self.spectrum[charge][key + 1]) == 0:
                        del self.spectrum[charge][key + 1]
                    get_idx_master_peak.add(get_master_peak_after)
                    get_idx_master_peak.add(peak)
                    self.merged += 1
                    if get_idx_master_peak.key() not in self.spectrum[charge]:
                        self.spectrum[charge][get_idx_master_peak.key()] = SortedList(key=lambda i: i.left)
                    self.spectrum[charge][get_idx_master_peak.key()].add(get_idx_master_peak)
        else:
            self.spectrum[charge][key] = SortedList(key=lambda i: i.left)
            self.add(peak, charge)

    def export_to_csv(self, path):
        with open(path, "w") as csvfile:
            writr = csv.writer(csvfile)  # , lineterminator=os.linesep)
            writr.writerow(
                (
                    "mz",
                    "intensity",
                    "counts",
                    "left border",
                    "right border",
                    "start_mz",
                    "ms1_charge",
                    "rel_intensity_ratio",
                    "counts_ratio",
                )
            )
            for charges in self.spectrum.keys():
                for key in self.spectrum[charges].keys():
                    for mp in self.spectrum[charges][key]:
                        writr.writerow(
                            (
                                mp.mz,
                                mp.intensity,
                                mp.counts,
                                mp.left,
                                mp.right,
                                mp.mz_origin,
                                charges,
                                mp.rel_intensity_ratio,
                                mp.counts_ratio,
                            )
                        )

    def load_from_tims(self, spectrum, ignore_charges: bool, delta_func=None):
        rel_int = _calculate_relative_intensity(spectrum["combined_INTENSITIES"])
        if delta_func is None:
            delta_func = _calculate_delta_by_ppm(40)
        for m, i in zip(spectrum["combined_MZ"], rel_int):
            p = Peak(float(m), float(i), delta_func)
            if ignore_charges:
                self.add(p, 0)
            else:
                raise NotImplementedError("Adding up intensities using precursor charge is not supported.")
