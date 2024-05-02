import pandas as pd
import logging

from pyopenms import IdXMLFile
from spectrum_fundamentals.mod_string import openms_to_internal, internal_without_mods
from tqdm import tqdm
from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)

def strToFloat(df):
  for col in df:
    try:
      df[col] = [float(i) for i in df[col]]
    except ValueError:
      continue
  return df

def readAndProcessIdXML(input_file, top=1):
  """
    Convert the (.idXML) format identification file to a DataFrame.

    :param input_file: Path to the input .idXML file.
    :param top: Number of top hits to consider, defaults to 1.
    :return: DataFrame containing identification information.
    """
  prot_ids = []; pep_ids = []
  IdXMLFile().load(input_file, prot_ids, pep_ids)
  meta_value_keys = []
  rows = []
  for peptide_id in pep_ids:
    spectrum_id = peptide_id.getMetaValue("spectrum_reference")
    scan_nr = spectrum_id[spectrum_id.rfind('=') + 1 : ] 

    hits = peptide_id.getHits()

    psm_index = 1
    for h in hits:
      if psm_index > top:
        break
      charge = h.getCharge() 
      score = h.getScore() 

      if "target" in h.getMetaValue("target_decoy"): 
          label = 1
      else:
          label = 0

      sequence = h.getSequence().toString() 

      if len(meta_value_keys) == 0: # fill meta value keys on first run
        h.getKeys(meta_value_keys)
        meta_value_keys = [x.decode() for x in meta_value_keys]

        all_columns = ['SpecId','PSMId','Label','Score','ScanNr','Peptide','peplen','ExpMass','charge','accessions'] + meta_value_keys

      # static part
      accessions = ';'.join([s.decode() for s in h.extractProteinAccessionsSet()])

      row = [spectrum_id, psm_index, label, score, scan_nr, sequence, str(len(sequence)), peptide_id.getMZ(), charge, accessions]
      # scores in meta values
      for k in meta_value_keys:
        s = h.getMetaValue(k)
        if type(s) == bytes:
          s = s.decode()
        row.append(s)
      rows.append(row)
      psm_index += 1
      break; # parse only first hit
  
  df =pd.DataFrame(rows, columns=all_columns)

  convert_dict = {'SpecId': str,
                  'PSMId': int,
                  'Score': float,
                  'ScanNr': int,
                  'peplen': int                
               }
  
  df = df.astype(convert_dict)
  df['Label'] = df['Label'].astype(bool)

  for prot_id in prot_ids:
     raw_file = (prot_id.getMetaValue("spectra_data"))[0].decode('utf-8').split('/')[-1].split('.')[0]
     break

  df['raw_file'] = raw_file 

  return df

class OpenMS(SearchResults):
    """Handle search results from OpenMS."""

    def read_result(self, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading OpenMS idXML file")

        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = list(self.path.rglob("*.idXML"))
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")
        
        openms_results = []
        for openms_file in tqdm(file_list):
            print("collecting data from ", openms_file)
            openms_results.append(readAndProcessIdXML(str(openms_file)))
        
        df = pd.concat(openms_results)

        logger.info("Finished reading OpenMS idXML file.")

        df = update_columns_for_prosit(df, tmt_labeled)
            
        return filter_valid_prosit_sequences(df)


def update_columns_for_prosit(df, tmt_labeled: str) -> pd.DataFrame:
    """
    Update columns of df to work with Prosit.

    :param df: df to modify
    :param tmt_labeled: True if tmt labeled
    :return: modified df as pd.DataFrame
    """
    if tmt_labeled != "":
        logger.debug("not implemented TMT modifications")
    else:
        df["MODIFIED_SEQUENCE"] = openms_to_internal(df["Peptide"].to_list())

    df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
    df["REVERSE"] = ~df["Label"]

    df.rename(
        columns={
            "raw_file": "RAW_FILE",
            "ExpMass": "MASS",
            "peplen": "PEPTIDE_LENGTH",
            "charge": "PRECURSOR_CHARGE",
            "ScanNr": "SCAN_NUMBER",
            "Score": "SCORE",
            "PSMId": "SCAN_EVENT_NUMBER" 
        },
        inplace=True,
    )

    # Select columns to return
    return_columns = ["RAW_FILE", "SCAN_NUMBER", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE", "SCAN_EVENT_NUMBER", "MASS", "SCORE", "REVERSE", "SEQUENCE", "PEPTIDE_LENGTH"]


    # if NA XL modification available
    if "NuXL:NA" in df.columns:
      df.rename(
          columns={
              "NuXL:NA": "NA_MOD",
          },
          inplace=True,
      )

      return_columns.append("NA_MOD")

    return df[return_columns]

