
extract_d = function(d_path, txt_path, path_to_bruker_dll, file_name){

    require(tidyverse)
    require(data.table)
    require(opentimsr)

    spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

    setup_bruker_so(path_to_bruker_dll)
    accept_Bruker_EULA_and_on_Windows_or_Linux <- TRUE

    if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
        path_to_bruker_dll <- path_to_bruker_dll
        setup_bruker_so(path_to_bruker_dll)
        all_columns <- c("frame", "scan", "tof", "intensity", "mz",
            "inv_ion_mobility", "retention_time")
    } else {
        all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
    }

    tbl_msms <- fread(paste(txt_path, "/msms.txt", sep = "")) %>%
        as_tibble() %>%
        rename_with(toupper) %>%
        spaceless() %>%
        filter(RAW_FILE == file_name)

    tbl_precursors <- fread(paste(txt_path,
        "/accumulatedMsmsScans.txt", sep = "")) %>%
        as_tibble() %>%
        rename_with(toupper) %>%
        spaceless() %>%
        filter(RAW_FILE == file_name) %>%
        filter(SCAN_NUMBER %in% tbl_msms$SCAN_NUMBER) %>%
        mutate(PRECURSOR = PASEF_PRECURSOR_IDS) %>%
        # Split the precursor ID and create a new row for each
        mutate(PRECURSOR = strsplit(as.character(PRECURSOR), ";")) %>%
        unnest(PRECURSOR) %>%
        distinct()

    scan_precursor_map <- tbl_precursors %>%
        select(SCAN_NUMBER, PRECURSOR) %>%
        distinct()

    tbl_pasef <- fread(paste(txt_path, "/pasefMsmsScans.txt", sep = "")) %>%
        as_tibble() %>%
        rename_with(toupper) %>%
        spaceless() %>%
        filter(RAW_FILE == file_name) %>%
        filter(PRECURSOR %in% tbl_precursors$PRECURSOR) %>%
        mutate(COLLISION_ENERGY = COLLISIONENERGY) %>%
        select(PRECURSOR, FRAME, SCANNUMBEGIN, SCANNUMEND, COLLISION_ENERGY) %>%
        distinct()

    df_tims <- OpenTIMS(d_path) # get data handle
    df_raw <- query(df_tims, tbl_pasef$FRAME, all_columns)

    tbl_raw <- df_raw %>%
        as_tibble() %>%
        rename_with(toupper) %>%
        spaceless() %>%
        distinct()

    dt_raw <- data.table(tbl_raw)
    dt_pasef <- data.table(tbl_pasef)

    dt_scan_precursor_map <- data.table(
        scan_precursor_map %>% mutate(PRECURSOR = as.integer(PRECURSOR)))

    dt_msms <- data.table(tbl_msms)
    return(list(dt_raw, dt_pasef, dt_scan_precursor_map, dt_msms))
  }