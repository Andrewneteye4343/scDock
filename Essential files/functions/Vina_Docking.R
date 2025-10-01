# functions/Vina_Docking.R
Vina_Docking <- function(Vina_Docking_input_path,
                         Vina_Docking_output_path,
                         Vina_Docking_ligand_ref_file = "/path/to/ligand_reference.csv",
                         Vina_Docking_receptor_ref_file = "path/to/receptor_reference.csv",
                         Vina_Docking_cas_txt_file = NULL,
                         Vina_Docking_use_fda = FALSE,
                         Vina_Docking_fda_txt = NULL,
                         Vina_Docking_docking_ligand_dir = NULL,
                         Vina_Docking_docking_receptor_dir = NULL,
                         Vina_Docking_vina_exhaustiveness = 8,
                         Vina_Docking_vina_num_modes = 9,
                         Vina_Docking_vina_seed = 42,
                         Vina_Docking_vina_cpu = 1
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  library(dplyr)

  # ---- grid.txt read function ----
  read_grid_file <- function(grid_file) {
    lines <- readLines(grid_file)
    vals <- sapply(lines, function(x) as.numeric(strsplit(x, "=")[[1]][2]))
    names(vals) <- sapply(lines, function(x) trimws(strsplit(x, "=")[[1]][1]))
    return(vals)
  }

  # ---- Vina docking running function ----
  run_vina_docking <- function(ligand_file, receptor_file, out_dir, params) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ligand_name <- tools::file_path_sans_ext(basename(ligand_file))
    output_pdbqt <- file.path(out_dir, paste0(ligand_name, "_AutoDockVina_result_structure.pdbqt"))
    log_file <- file.path(out_dir, paste0(ligand_name, "_AutoDockVina_result_score.txt"))

    grid_file <- sub("\\.pdbqt$", "_grid.txt", receptor_file)
    if (!file.exists(grid_file)) {
      stop("[Vina] Missing grid file: ", grid_file, " for receptor: ", receptor_file)
    }
    grid <- read_grid_file(grid_file)

    cmd <- c("--ligand", ligand_file,
             "--receptor", receptor_file,
             "--out", output_pdbqt,
             "--log", log_file,
             "--exhaustiveness", as.character(params$exhaustiveness),
             "--num_modes", as.character(params$num_modes),
             "--seed", as.character(params$seed),
             "--cpu", as.character(params$cpu),
             "--center_x", as.character(grid["center_x"]),
             "--center_y", as.character(grid["center_y"]),
             "--center_z", as.character(grid["center_z"]),
             "--size_x", as.character(grid["size_x"]),
             "--size_y", as.character(grid["size_y"]),
             "--size_z", as.character(grid["size_z"])
    )
    message("")
    message("[Vina] Docking ligand ", ligand_name, " with receptor ", basename(receptor_file))
    status <- system2("vina", args = cmd)
    if (status != 0) {
      warning("[Vina] Docking failed for ligand ", ligand_file, " with receptor ", receptor_file)
    }
  }

  # ---- Set Vina arguments ----
  vina_params <- list(
    exhaustiveness = Vina_Docking_vina_exhaustiveness,
    num_modes = Vina_Docking_vina_num_modes,
    seed = Vina_Docking_vina_seed,
    cpu = Vina_Docking_vina_cpu
  )

  # ---- Read LR CSV ----
  csv_files <- list.files(Vina_Docking_input_path, pattern = "^top\\d+_LR_.*\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    # fallback: try multi-group file
    multi_csv <- file.path(Vina_Docking_input_path, "multiGroup_significant_LR_by_prob_diff.csv")
    if (file.exists(multi_csv)) {
      message("[INFO] No top*_LR_*.csv found. Using multiGroup_significant_LR_by_prob_diff.csv instead.")
      csv_files <- c(multi_csv)
    } else {
      stop("No LR CSV files found (neither top*_LR_*.csv nor multiGroup_significant_LR_by_prob_diff.csv) in: ",
           Vina_Docking_input_path)
    }
  }

  ligands <- receptors <- c()
  for (csv in csv_files) {
    df <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next

    # Ensure ligand/receptor columns exist
    if (!all(c("ligand","receptor") %in% colnames(df)) && "interaction" %in% colnames(df)) {
      parts <- strsplit(df$interaction, "_")
      df$ligand <- sapply(parts, `[`, 1)
      df$receptor <- sapply(parts, `[`, 2)
    }

    if (all(c("ligand","receptor") %in% colnames(df))) {
      ligands <- c(ligands, df$ligand)
      receptors <- c(receptors, df$receptor)
    } else {
      warning("[Vina] CSV ", basename(csv), " does not contain ligand/receptor information.")
    }
  }
  ligands <- unique(ligands)
  receptors <- unique(receptors)

  # ---- Load reference ----
  ligand_ref <- read.csv(Vina_Docking_ligand_ref_file, stringsAsFactors = FALSE)
  receptor_ref <- read.csv(Vina_Docking_receptor_ref_file, stringsAsFactors = FALSE)
  ligand_match <- ligand_ref %>% filter(protein_name %in% ligands)
  receptor_match <- receptor_ref %>% filter(protein_name %in% receptors)

  if (!dir.exists(Vina_Docking_output_path)) dir.create(Vina_Docking_output_path, recursive = TRUE)
  write.csv(ligand_match, file.path(Vina_Docking_output_path, "ligands_with_PDB.csv"), row.names = FALSE)
  write.csv(receptor_match, file.path(Vina_Docking_output_path, "receptors_with_PDB.csv"), row.names = FALSE)

  # ---- Create PDB directory ----
  ligand_dir <- file.path(Vina_Docking_output_path, "ligand_from_PDB_LR_pairs")
  receptor_dir <- file.path(Vina_Docking_output_path, "receptor_from_PDB_LR_pairs")
  if (!dir.exists(ligand_dir)) dir.create(ligand_dir, recursive = TRUE)
  if (!dir.exists(receptor_dir)) dir.create(receptor_dir, recursive = TRUE)

  # ---- Download protein structure ----
  download_structures <- function(df, outdir, label) {
    for (raw_ids in df$PDB_model) {
      if (is.na(raw_ids) || raw_ids=="") next
      ids <- trimws(unlist(strsplit(raw_ids, ";")))
      for (id in ids) {
        if (id=="") next
        subdir <- file.path(outdir, id)
        if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
        outfile <- file.path(subdir, paste0(id, ".pdb"))

        if (grepl("^AF-", id)) {
          if (!file.exists(outfile)) {
            message("[",label,"] Downloading AlphaFold model: ", id)
            system2("python3", args = c(normalizePath("functions/download_alphafold.py"), id, subdir))
          }
          message("[",label,"] Preprocessing with AutoDockTools: ", outfile)
          system2("python3", args = c(normalizePath("functions/prepare_receptor.py"), outfile, subdir))
        } else if (grepl("^[0-9][A-Za-z0-9]{3}\\.[A-Z]$", id)) {
          parts <- strsplit(id, "\\.")[[1]]
          pdb_id <- tolower(parts[1]); chain <- toupper(parts[2])
          outfile <- file.path(subdir, paste0(pdb_id,".1.",chain,".pdb"))
          if (!file.exists(outfile)) {
            message("[",label,"] Downloading SWISS-MODEL: ", pdb_id, " chain ", chain)
            system2("python3", args = c(normalizePath("functions/download_swissmodel.py"), pdb_id, chain, subdir))
          }
          message("[",label,"] Preprocessing with AutoDockTools: ", outfile)
          system2("python3", args = c(normalizePath("functions/prepare_receptor.py"), outfile, subdir))
        } else warning("[",label,"] Unrecognized ID format: ", id)
      }
    }
  }

  download_structures(ligand_match, ligand_dir, "ligand")

  # only download receptors if user did NOT provide Vina_Docking_docking_receptor_dir
  if (is.null(Vina_Docking_docking_receptor_dir)) {
    download_structures(receptor_match, receptor_dir, "receptor")
  }

  # ---- Ligand structure source ----
  cas_dir <- NULL

  if (Vina_Docking_use_fda) {
    # ---- FDA compounds support ----
    if (is.null(Vina_Docking_fda_txt) || !file.exists(Vina_Docking_fda_txt)) {
      stop("[FDA] Vina_Docking_use_fda = TRUE but Vina_Docking_fda_txt is missing or does not exist!")
    }

    fda_dir <- file.path(Vina_Docking_output_path, "fda_results")
    if (!dir.exists(fda_dir)) dir.create(fda_dir, recursive = TRUE)

    fda_cas <- trimws(readLines(Vina_Docking_fda_txt))
    if (length(fda_cas) > 0) {
      message("[FDA] Downloading structures for ", length(fda_cas), " FDA compounds...")
      system2("python3", args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        Vina_Docking_fda_txt, fda_dir
      ))
      Vina_Docking_docking_ligand_dir <- fda_dir
    } else {
      stop("[FDA] fda.txt is empty! Please provide CAS numbers.")
    }

  } else {
    # ---- CAS compounds support ----
    if (!is.null(Vina_Docking_cas_txt_file) && file.exists(Vina_Docking_cas_txt_file)) {
      cas_dir <- file.path(Vina_Docking_output_path, "ligand_structures_from_CAStxt_for_AutoDockVina")
      if (!dir.exists(cas_dir)) dir.create(cas_dir)
      cas_numbers <- trimws(readLines(Vina_Docking_cas_txt_file))
      message("[CAS] Downloading structures for ", length(cas_numbers), " compounds...")
      system2("python3", args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        Vina_Docking_cas_txt_file, cas_dir
      ))
      Vina_Docking_docking_ligand_dir <- cas_dir
    }
  }

  # ---- AutoDock Vina docking ----
  Vina_Docking_docking_ligand_dir <- if (is.null(Vina_Docking_docking_ligand_dir)) cas_dir else Vina_Docking_docking_ligand_dir

  if (is.null(Vina_Docking_docking_receptor_dir)) {
    Vina_Docking_docking_receptor_dir <- list(ligand_dir, receptor_dir)
  } else {
    message("[INFO] Using user-provided receptor directory: ", Vina_Docking_docking_receptor_dir)
  }

  receptor_files <- unlist(lapply(Vina_Docking_docking_receptor_dir, function(d) {
    list.files(d, pattern = "_prepared\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  }))

  ligand_files <- list.files(Vina_Docking_docking_ligand_dir, pattern = "\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  ligand_files <- ligand_files[file.size(ligand_files) > 0]

  for (receptor in receptor_files) {
    rec_subdir <- dirname(receptor)
    for (ligand in ligand_files) {
      run_vina_docking(
        ligand_file = ligand,
        receptor_file = receptor,
        out_dir = rec_subdir,
        params = vina_params
      )
    }

    # ---- Collect docking results per receptor ----
    message("[Vina] Collecting docking results for receptor folder: ", rec_subdir)
    log_files <- list.files(rec_subdir,
                            pattern = "_AutoDockVina_result_score.txt$",
                            full.names = TRUE,
                            recursive = FALSE)

    results <- lapply(log_files, function(logf) {
      lines <- readLines(logf, warn = FALSE)
      affinity <- NA
      # Get the first mode 
      line1 <- grep("^\\s*1\\s", lines, value = TRUE)
      if (length(line1) > 0) {
        affinity <- as.numeric(strsplit(trimws(line1), "\\s+")[[1]][2])
      }
      ligand <- sub("_AutoDockVina_result_score.txt$", "", basename(logf))
      data.frame(ligand = ligand,
                 affinity = affinity,
                 stringsAsFactors = FALSE)
    })

    if (length(results) > 0) {
      results_df <- do.call(rbind, results) %>%
        arrange(affinity)
      out_csv <- file.path(rec_subdir, "AutoDockVina_score.csv")
      write.csv(results_df, out_csv, row.names = FALSE)
      message("[Vina] Receptor score summary saved to: ", out_csv)
    }
  }
}