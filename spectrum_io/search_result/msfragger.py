import logging
import re

import pandas as pd
import spectrum_fundamentals.constants as c

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    @staticmethod
    def read_result(path: str, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param path: path to msms.txt to read
        :param tmt_labeled: tmt label as str
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading msfragger tsv file")
        df = pd.read_csv(
            path,
            usecols=lambda x: x.upper()
            in [
                "PEPTIDE",
                "PROTEIN",
                "PEPTIDE LENGTH",
                "SPECTRUM FILE",
                "SPECTRUM",
                "ASSIGNED MODIFICATIONS",
                "CHARGE",
                "LABELING STATE",
                "OBSERVED MASS",  # = Calculated Precursor mass; TODO get column with experimental Precursor mass instead
                "HYPERSCORE",
                "RETENTION",
            ],
            sep="\t",
        )
        logger.info("Finished reading msfragger tsv file")

        df.rename(
            columns={
                "Peptide": "SEQUENCE",
                "Assigned Modifications": "MODIFICATIONS",
                "Observed Mass": "MASS",
                "Hyperscore": "SCORE",
                "Retention": "RETENTION TIME",
                "Spectrum File": "RAW FILE",
            },
            inplace=True,
        )

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")
        df.rename(columns={"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
        df[["RAW_FILE", "SCAN_NUMBER"]] = df["SPECTRUM"].str.split(".", expand=True, n=2)[[0, 1]]
        df["REVERSE"] = df["PROTEIN"].str.contains("Reverse")

        df["MODIFICATIONS"] = df["MODIFICATIONS"].fillna(0)
        mod_masses_reverse = {round(float(v), 3): k for k, v in c.MOD_MASSES.items()}
        sequences = []
        for _, row in df.iterrows():
            modifications = row["MODIFICATIONS"]
            if modifications == 0:
                sequences.append(row["SEQUENCE"])
            else:
                modifications = modifications.split(", ")
                sequence = row["SEQUENCE"]
                skip = 0
                for mod in sorted(
                    modifications, key=lambda s: 0 if s.startswith("N") else int(re.sub("[^0-9]", "", s.split("(")[0]))
                ):
                    pos, mass = mod.split("(")
                    mass = mass.replace(")", "")
                    if pos == "N-term":
                        continue
                    else:
                        pos = re.sub("[^0-9]", "", pos)
                    sequence = (
                        sequence[: int(pos) + skip]
                        + mod_masses_reverse[round(float(mass), 3)]
                        + sequence[int(pos) + skip :]
                    )
                    skip = skip + len(mod_masses_reverse[round(float(mass), 3)])
                sequences.append(sequence)
        df["MODIFIED_SEQUENCE"] = sequences

        logger.info(f"No of sequences before Filtering is {len(df['PEPTIDE_LENGTH'])}")
        df = df[(df["PEPTIDE_LENGTH"] <= 30)]
        df = df[(~df["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))]
        df = df[(~df["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))]
        df = df[(~df["SEQUENCE"].str.contains("U"))]
        df = df[df["PRECURSOR_CHARGE"] <= 6]
        df = df[df["PEPTIDE_LENGTH"] >= 7]
        logger.info(f"No of sequences after Filtering is {len(df['PEPTIDE_LENGTH'])}")
        return df
