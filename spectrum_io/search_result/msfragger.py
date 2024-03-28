import logging
from pathlib import Path
from typing import Union

import pandas as pd
import spectrum_fundamentals.constants as c
from pyteomics import pepxml
from spectrum_fundamentals.mod_string import internal_without_mods, msfragger_to_internal
from tqdm import tqdm

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    def read_result(self, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :return: pd.DataFrame with the formatted data
        """
        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = list(self.path.rglob("*.pepXML"))
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")

        mod_table = pd.DataFrame()
    	#just in case  not in all pepxml files are the same modification lines at the start of the file
        for pepxml_file in file_list:
        	new_mod_table = read_out_pepxml_modification_table(pepxml_file)
        	mod_table = pd.concat([mod_table, new_mod_table], ignore_index=True)
        	mod_table = mod_table.drop_duplicates()
        
        ms_frag_results = []
        for pep_xml_file in tqdm(file_list):
            ms_frag_results.append(pepxml.DataFrame(str(pep_xml_file)))

        df = pd.concat(ms_frag_results)

        df = add_fpop_and_existing_status_to_df(df=df, unimod_table=mod_table)
        
        df = df[df["is_existing_mod"]==True]

        df = update_columns_for_prosit(df, tmt_labeled)

        return df
        #return filter_valid_prosit_sequences(df)


def update_columns_for_prosit(df, tmt_labeled: str) -> pd.DataFrame:
    """
    Update columns of df to work with Prosit.

    :param df: df to modify
    :param tmt_labeled: True if tmt labeled
    :return: modified df as pd.DataFrame
    """
    df["REVERSE"] = df["protein"].apply(lambda x: "rev" in str(x))
    df["RAW_FILE"] = df["spectrum"].apply(lambda x: x.split(".")[0])
    df["MASS"] = df["precursor_neutral_mass"]
    df["PEPTIDE_LENGTH"] = df["peptide"].apply(lambda x: len(x))

    if tmt_labeled != "":
        unimod_tag = c.TMT_MODS[tmt_labeled]
        logger.info("Adding TMT fixed modifications")
        df["MODIFIED_SEQUENCE"] = msfragger_to_internal(
            df["modified_peptide"].to_list(),
            fixed_mods={"C": "C[UNIMOD:4]", r"n[\d+]": f"{unimod_tag}-", "K": f"K{unimod_tag}"},
        )
    else:
        df["MODIFIED_SEQUENCE"] = msfragger_to_internal(df["modified_peptide"].to_list())

    df.rename(
        columns={
            "assumed_charge": "PRECURSOR_CHARGE",
            "index": "SCAN_EVENT_NUMBER",
            "peptide": "SEQUENCE",
            "start_scan": "SCAN_NUMBER",
            "hyperscore": "SCORE",
        },
        inplace=True,
    )
    df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
    return df[
        [
            "RAW_FILE",
            "SCAN_NUMBER",
            "MODIFIED_SEQUENCE",
            "PRECURSOR_CHARGE",
            "SCAN_EVENT_NUMBER",
            "MASS",
            "SCORE",
            "REVERSE",
            "SEQUENCE",
            "PEPTIDE_LENGTH",
        ]
    ]

def read_out_pepxml_modification_table(pepxml_file):
    ####! with changed molecular weights to account for slight differences in the number of digits
    """
    Creates a table with the different possible aminoacid modifications
    based on the aa modification lines at the start of the .pepxml file
    Also columns for the different modification names like C[160] (for modifications displayed in pepxml) 
    or C[57] (for modifications displayed in psm.tsv) are added
    """
    # Read in amino modifications from .pepxml
    modification_lines = []
    with open(pepxml_file, "r") as file:
        for line in file:
            if line.startswith("<aminoacid_modification"):
                modification_lines.append(line.strip())


    # Create modification dataframe with massdifferences and mass
    aminoacids = []
    mass_values = []
    mass_diff_values = []

    pattern_aa = r'aminoacid="(\w)"'
    pattern_mass = r'mass="([-]?[\d.]+)"'
    pattern_mass_diff = r'massdiff="([-]?[\d.]+)"'

    # Iterate through the input strings and extract data
    for input_string in modification_lines:
        match = re.search(pattern_aa, input_string)
        aminoacid_letter = match.group(1)
        aminoacids.append(aminoacid_letter)

        match = re.search(pattern_mass_diff, input_string)
        mass_diff_value = float(match.group(1))
        mass_diff_values.append(mass_diff_value)

        match = re.search(pattern_mass, input_string)
        mass_value = float(match.group(1))
        mass_values.append(mass_value)


    df_pep_xml_aa_mods = pd.DataFrame({"Aminoacid": aminoacids, "Massdifference": mass_diff_values, "Mass": mass_values})

    # Add looked up / identified unimod modifications
    # in case of key errors  because of the mass differences, like for 57.02146(1), add additional entries
    unimod_numbers={"57.021461": "UNIMOD:4",
                    "57.02146": "UNIMOD:4",
                    "15.9949": "UNIMOD:35",
                    "31.989829": "UNIMOD:425",
                    "31.989828": "UNIMOD:425",
                    "47.984744": "UNIMOD:345",
                    "47.984745": "UNIMOD:345",
                    "13.979265": "UNIMOD:1918",
                    "-43.053433": "UNIMOD:344",
                    "-43.053432": "UNIMOD:344",
                    "-22.031969": "UNIMOD:349",
                    "-22.03197": "UNIMOD:349",
                    "-23.015984": "UNIMOD:348",
                    "-10.031969": "UNIMOD:1916",
                    "4.9735": "UNIMOD:1917",
                    "-30.010565": "UNIMOD:1915",
                    "-27.994915": "UNIMOD:369",
                    "-43.989829": "UNIMOD:553",
                    "-43.98983": "UNIMOD:553",
                    "-25.031631": "UNIMOD:10000",
                    "-9.036716": "UNIMOD:10001"}
    def lookup_unimod_name(row):
        aa = row["Aminoacid"]
        massdiff_value = str(row["Massdifference"])
        unimod_name = aa+"["+unimod_numbers[massdiff_value]+"]"
        return unimod_name
    
    df_pep_xml_aa_mods["Unimod_name"] = df_pep_xml_aa_mods.apply(lookup_unimod_name, axis=1)
    
    # round mass values
    df_pep_xml_aa_mods_round = df_pep_xml_aa_mods.copy()
    df_pep_xml_aa_mods_round["Mass"] = df_pep_xml_aa_mods_round["Mass"].round().astype(int)
    df_pep_xml_aa_mods_round["Massdifference"] = df_pep_xml_aa_mods_round["Massdifference"].round().astype(int)

    # Add looked up / identified modifications
    modifications_from_unimod={"57": "Iodoacetamide derivative",
                           "16": "Oxidation",
                           "32": "Dioxidation",
                           "48": "Trioxidation",
                           "14": "Carbonylation",
                           "-43": "Arg deguanidation",
                           "-22": "H to D",
                           "-23": "H to N",
                           "-10": "H - 10",
                           "5": "H + 5",
                           "-30": "Decarboxylation",
                           "-28": "CO loss",
                           "-44": "CO2 loss",
                           "-25": "C dioxidation assuming fixed carbamidomethylation",
                           "-9": "C trioxidation assuming fixed carbamidomethylation"}
    def lookup_mod(row):
        massdiff_value = str(row["Massdifference"])
        modification = modifications_from_unimod[massdiff_value]
        return modification

    df_pep_xml_aa_mods_round["Modification"] = df_pep_xml_aa_mods_round.apply(lookup_mod, axis=1)

    #adding apepxml&psm column (like "C[160]")
    def create_pepxml_mod_name(row):
        aa = row["Aminoacid"]
        mass = str(row["Mass"])
        pepxml_mod_name = aa+"["+mass+"]"
        return pepxml_mod_name

    def create_psm_mod_name(row):
        aa = row["Aminoacid"]
        mass = str(row["Massdifference"])
        psm_mod_name = aa+"["+mass+"]"
        return psm_mod_name

    df_pep_xml_aa_mods_round["pepxml_mod_name"] = df_pep_xml_aa_mods_round.apply(create_pepxml_mod_name, axis=1)
    df_pep_xml_aa_mods_round["psm_mod_name"] = df_pep_xml_aa_mods_round.apply(create_psm_mod_name, axis=1)


    # Add boolean value if modification exists according to msfragger paper (table 2)
    exisiting_mods = {"57": "C",
                  "16": "MFHILVWYADEKNPQR",
                  "32": "CFMWY",
                  "48": "CFWY",
                  "14": "EIKLPQRV",
                  "-43": "R",
                  "-22": "H",
                  "-23": "H",
                  "-10": "H",
                  "5": "H",
                  "-30": "DE",
                  "-28": "DE",
                  "-44": "DE",
                  "-25": "C",
                  "-9":"C"
    }
    def add_exisiting_value(row):
        massdiff = str(row["Massdifference"])
        aa = row["Aminoacid"]
        if massdiff in exisiting_mods.keys():
            bool_val = aa in exisiting_mods[massdiff]
        else:
            bool_val = False
        return bool_val
    df_pep_xml_aa_mods_round["is_existing_mod"] = df_pep_xml_aa_mods_round.apply(add_exisiting_value, axis=1)
    
    def add_fpop_status_unimod(row):
        list_of_oxidation_mod = ["Oxidation", "Dioxidation", "Trioxidation",
                                    "C dioxidation assuming fixed carbamidomethylation", "C trioxidation assuming fixed carbamidomethylation" ]
        if row["Modification"] in list_of_oxidation_mod:
            return True
        else:
            return False
    df_pep_xml_aa_mods_round["fpop_mod"] = df_pep_xml_aa_mods_round.apply(add_fpop_status_unimod, axis=1)

    return df_pep_xml_aa_mods_round
    
def add_fpop_and_existing_status_to_df(df, unimod_table=pd.DataFrame(), asign_modus: Union["one_non_exist_to_exclude", "one_exist_to_include"] = "one_non_exist_to_exclude", keep_unmodified=True):
    """
    For adding a "fpop_mod" and "is_existing_mod" column to a pyteomics.pepxml dataframe.
    "fpop_mod" is set True/1 if at least one of the modifications in the PSM are counted as a oxidation modification.
    The setting "is_existing_mod" is dependend of the asign_modus:
    If "one_non_exist_to_exclude" is set as asign-modus, than for every PSM, that has at least one modification, that shouldn't exist,
    according to the fpop paper Table 2, the "is_existing_mod" is set as False.
    If "one_exist_to_include" is set as asign-modus, than for every PSM, that has at least one modification, that should exist,
    according to the fpop paper Table 2, the "is_existing_mod" is set as True.
    """
    new_df = df.copy()
    new_df["fpop_mod"] = False
    new_df["is_existing_mod"] = False
    if "modified_peptide" in new_df.columns:
        modus = "pepxml"
    if "Assigned Modifications" in new_df.columns:
        modus = "psm"
    if unimod_table.empty:
        unimod_table = create_modification_table()
    unimod_table = unimod_table.set_index(modus+"_mod_name")
    new_df["fpop_mod"] = new_df.apply(lambda x: add_fpop_status_to_row(x,modus=modus,unimod_table=unimod_table), axis=1)
    new_df["is_existing_mod"] = new_df.apply(lambda x: add_existing_status_to_row(x,modus=modus,unimod_table=unimod_table, asign_modus=asign_modus, keep_unmodified=keep_unmodified), axis=1)
    return new_df

def add_fpop_status_to_row(row, modus: Union["pepxml", "psm"], unimod_table):
    if modus == "pepxml":
        peptide = row["modified_peptide"]
        possible_mods = r"\w\[\d+\]"
        matches = re.findall(possible_mods, peptide)
        if matches:
            for match in matches:
                if unimod_table.loc[match]["fpop_mod"]:
                    return unimod_table.loc[match]["fpop_mod"] 
            return False
        else:
            return False
    if modus == "psm":
        modifications_row = row["Assigned Modifications"]
        if not (isinstance(modifications_row, (float, int)) and math.isnan(modifications_row)):
            possible_mods = r'[A-Za-z]\(-?\d+\.\d+\)'
            matches = re.findall(possible_mods, modifications_row)
            for match in matches:
                massdiff = round(float(match[2:-1]))
                aa = match[0]
                mod = aa+"["+str(massdiff)+"]"
                if unimod_table.loc[mod]["fpop_mod"]:
                    return unimod_table.loc[mod]["fpop_mod"]
            return False
        return False
    
def add_existing_status_to_row(row, modus: Union["pepxml", "psm"], unimod_table, asign_modus: Union["one_non_exist_to_exclude", "one_exist_to_include"], keep_unmodified=True):
    if modus == "pepxml":
        peptide = row["modified_peptide"]
        possible_mods = r"\w\[\d+\]"
        matches = re.findall(possible_mods, peptide)
        if matches:
            for match in matches:
                if asign_modus == "one_exist_to_include":
                    if unimod_table.loc[match]["is_existing_mod"]:
                        # in this case the entry in the table is already True
                        return  True
                if asign_modus == "one_non_exist_to_exclude":
                    if not unimod_table.loc[match]["is_existing_mod"]:
                        # in this case the entry in the table is already False
                        return False
            if asign_modus == "one_exist_to_include":
                return False
            if asign_modus == "one_non_exist_to_exclude":
                return True
        else:
            if keep_unmodified:
                return True
            else:
                return False
    if modus == "psm":
        modifications_row = row["Assigned Modifications"]
        # to account for the possible NaN values. Without the isinstance() ckeck, error "must be real number, not str" can occur
        if not (isinstance(modifications_row, (float, int)) and math.isnan(modifications_row)):
            possible_mods = r'[A-Za-z]\(-?\d+\.\d+\)'
            matches = re.findall(possible_mods, modifications_row)
            for match in matches:
                massdiff = round(float(match[2:-1]))
                aa = match[0]
                mod = aa+"["+str(massdiff)+"]"
                if asign_modus == "one_exist_to_include":
                    if unimod_table.loc[mod]["is_existing_mod"]:
                        return True
                if asign_modus == "one_non_exist_to_exclude":
                    if not unimod_table.loc[mod]["is_existing_mod"]:
                        return False
            if asign_modus == "one_exist_to_include":
                return False
            if asign_modus == "one_non_exist_to_exclude":
                return True
        else: # = no annotated modifications 
            if keep_unmodified:
                return True
            else:
                return False
