from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from pyopenms import IdXMLFile, PeptideIdentification, ProteinIdentification
from spectrum_fundamentals.mod_string import internal_without_mods
from tqdm import tqdm

from .search_results import SearchResults, parse_mods

logger = logging.getLogger(__name__)


def _extract_scan_number(spectrum_id: int | float | bytes | str | list[int] | list[float] | list[bytes]) -> int:
    """Extract scan number from spectrum ID."""
    if isinstance(spectrum_id, str):
        return int(spectrum_id[spectrum_id.rfind("=") + 1 :])
    elif isinstance(spectrum_id, bytes):
        decoded = spectrum_id.decode("utf-8")
        return int(decoded[decoded.rfind("=") + 1 :])
    else:
        raise TypeError(f"spectrum_reference must be a str or bytes, but got {type(spectrum_id)}")


def _get_raw_file_name(prot_ids: list[ProteinIdentification]) -> str:
    """Extract raw file name from the first ProteinIdentification."""
    for prot_id in prot_ids:
        spectra_data = prot_id.getMetaValue("spectra_data")
        if isinstance(spectra_data, list) and spectra_data:
            first_element = spectra_data[0]
            if isinstance(first_element, bytes):
                return first_element.decode("utf-8").split("/")[-1].split(".")[0]
            else:
                raise TypeError(f"Expected bytes in the list, but got {type(first_element)}")
        else:
            raise TypeError("spectra_data must be a non-empty list of bytes")
    raise ValueError("Raw file name could not be extracted from protein identifications.")


def _read_and_process_id_xml(input_file: Path, top: int = 0) -> pd.DataFrame:
    """
    Convert the (.idXML) format identification file to a DataFrame.

    :param input_file: Path to the input .idXML file.
    :param top: Number of top hits to consider, defaults to 0, which returns all hits.
    :return: DataFrame containing identification information.
    """
    prot_ids: list[ProteinIdentification] = []
    pep_ids: list[PeptideIdentification] = []
    IdXMLFile().load(str(input_file), prot_ids, pep_ids)

    meta_value_keys_bytes: list[bytes] = []
    meta_value_keys: list[str] = []
    rows = []
    for peptide_id in pep_ids:

        spectrum_id = peptide_id.getMetaValue("spectrum_reference")
        # extract scan number
        scan_nr = _extract_scan_number(spectrum_id=spectrum_id)

        hits = peptide_id.getHits()
        psm_index = 1
        for h in hits:
            if top > 0 and psm_index > top:
                break
            charge = h.getCharge()
            score = h.getScore()

            if "target" in str(h.getMetaValue("target_decoy")):
                label = 1
            else:
                label = 0

            sequence = h.getSequence().toString()

            if len(meta_value_keys_bytes) == 0:  # fill meta value keys on first run
                h.getKeys(meta_value_keys_bytes)
                meta_value_keys = [x.decode() for x in meta_value_keys_bytes]

                all_columns = [
                    "SpecId",
                    "PSMId",
                    "Label",
                    "Score",
                    "ScanNr",
                    "Peptide",
                    "peplen",
                    "ExpMass",
                    "charge",
                    "accessions",
                ] + meta_value_keys

            # static part
            accessions = ";".join(sorted([s.decode() for s in h.extractProteinAccessionsSet()]))

            row = [
                spectrum_id,
                psm_index,
                label,
                score,
                scan_nr,
                sequence,
                str(len(sequence)),
                peptide_id.getMZ(),
                charge,
                accessions,
            ]

            # scores in meta values
            for k in meta_value_keys:
                s = h.getMetaValue(k)
                if isinstance(s, bytes):
                    s = s.decode()
                row.append(s)
            rows.append(row)
            psm_index += 1

    df = pd.DataFrame(rows, columns=all_columns)

    df = df.astype(
        {
            "SpecId": str,
            "PSMId": "int64",
            "Score": float,
            "ScanNr": "int64",
            "peplen": "int64",
            "Label": bool,
            "charge": "int64",
        }
    )

    # extract raw file name
    raw_file = _get_raw_file_name(prot_ids)
    df["raw_file"] = raw_file

    return df


class OpenMS(SearchResults):
    """Handle search results from OpenMS."""

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {"C(Carbamidomethyl)": 4, "M(Oxidation)": 35, "R(Deamidated)": 7, "Q(Deamidated)": 7, "N(Deamidated)": 7}

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_label: tmt label as str
        :param custom_mods: dict with custom variable and static identifier and respective internal equivalent and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :raises NotImplementedError: in case TMT or ptm_unimod_id/ptm_sites are given.

        :return: pd.DataFrame with the formatted data
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))
        if tmt_label:
            raise NotImplementedError("TMT data is currently not supported for OpenMS")
            # unimod_tag = c.TMT_MODS[tmt_label]
            # parsed_mods[r"K\[\+\d+\.\d+\]"] = f"K{unimod_tag}"
            # parsed_mods[r"^\[\+\d+\.\d+\]"] = f"{unimod_tag}"

        if ptm_unimod_id is not None and ptm_unimod_id != 0 or ptm_sites is not None:
            raise NotImplementedError("ptm_unimod_id and ptm_sties are not supported yet.")

        logger.info("Reading OpenMS idXML file(s)...")

        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = list(self.path.rglob("*.idXML"))
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")

        openms_results = []
        for openms_file in tqdm(file_list):
            print("Collecting data from ", openms_file)
            openms_results.append(_read_and_process_id_xml(openms_file))

        self.results = pd.concat(openms_results)

        logger.info("Finished reading OpenMS idXML file(s).")

        self.convert_to_internal(mods=parsed_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites)

        return self.filter_valid_prosit_sequences()

    def filter_valid_prosit_sequences(self):
        """Filter valid Prosit sequences."""
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results.index)}")
        # retain only peptides that fall within [7, 30] length supported by Prosit
        self.results = self.results[(self.results["PEPTIDE_LENGTH"] <= 30) & (self.results["PEPTIDE_LENGTH"] >= 7)]
        # remove unsupported mods to exclude
        self.results = self.results[~self.results["MODIFIED_SEQUENCE"].str.contains(r"\[\d+\]", regex=True)]
        # remove precursor charges greater than 6
        self.results = self.results[self.results["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(self.results.index)}")

        return self.results

    def convert_to_internal(self, mods: dict[str, str], ptm_unimod_id: int | None, ptm_sites: list[str] | None):
        """
        Convert all columns in the Sage output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping Sage-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        df = self.results
        df.replace({"Peptide": mods}, regex=True, inplace=True)
        df["SEQUENCE"] = internal_without_mods(df["Peptide"])
        df["Label"] = ~df["Label"]

        df.rename(
            columns={
                "raw_file": "RAW_FILE",
                "ExpMass": "MASS",
                "peplen": "PEPTIDE_LENGTH",
                "charge": "PRECURSOR_CHARGE",
                "ScanNr": "SCAN_NUMBER",
                "Peptide": "MODIFIED_SEQUENCE",
                "Score": "SCORE",
                "PSMId": "SCAN_EVENT_NUMBER",
                "accessions": "PROTEINS",
                "Label": "REVERSE",
            },
            inplace=True,
        )

        # if NA XL modification available
        if "NuXL:NA" in df.columns:
            df.rename(
                columns={
                    "NuXL:NA": "NA_MOD",
                },
                inplace=True,
            )
