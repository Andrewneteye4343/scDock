vina_docking <- function(input_dir,
                         ligand_ref_file = "ligand_reference.csv",
                         receptor_ref_file = "receptor_reference.csv",
                         output_dir,
                         python_bin = "python3",
                         cas_txt_file = NULL,
                         docking_ligand_dir = NULL,
                         use_fda = FALSE,
                         fda_txt_path = NULL,
                         docking_receptor_dir = NULL,
                         vina_exhaustiveness = 8,
                         vina_num_modes = 9,
                         vina_seed = 42,
                         vina_cpu = 1
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
    output_pdbqt <- file.path(out_dir, paste0(ligand_name, "_out.pdbqt"))
    log_file <- file.path(out_dir, paste0(ligand_name, "_log.txt"))

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

    message("[Vina] Docking ligand ", ligand_name, " with receptor ", basename(receptor_file))
    status <- system2("vina", args = cmd)
    if (status != 0) {
      warning("[Vina] Docking failed for ligand ", ligand_file, " with receptor ", receptor_file)
    }
  }

  # ---- Set Vina arguments ----
  vina_params <- list(
    exhaustiveness = vina_exhaustiveness,
    num_modes = vina_num_modes,
    seed = vina_seed,
    cpu = vina_cpu
  )

  # ---- Read LR CSV ----
  csv_files <- list.files(input_dir, pattern = "^top\\d+_LR_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No LR CSV files found in: ", input_dir)

  ligands <- receptors <- c()
  for (csv in csv_files) {
    df <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && all(c("ligand","receptor") %in% colnames(df))) {
      ligands <- c(ligands, df$ligand)
      receptors <- c(receptors, df$receptor)
    }
  }
  ligands <- unique(ligands)
  receptors <- unique(receptors)

  # ---- Load reference ----
  ligand_ref <- read.csv(ligand_ref_file, stringsAsFactors = FALSE)
  receptor_ref <- read.csv(receptor_ref_file, stringsAsFactors = FALSE)
  ligand_match <- ligand_ref %>% filter(protein_name %in% ligands)
  receptor_match <- receptor_ref %>% filter(protein_name %in% receptors)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(ligand_match, file.path(output_dir, "ligands_with_PDB.csv"), row.names = FALSE)
  write.csv(receptor_match, file.path(output_dir, "receptors_with_PDB.csv"), row.names = FALSE)

  # ---- Create PDB directory ----
  ligand_dir <- file.path(output_dir, "ligand_PDB")
  receptor_dir <- file.path(output_dir, "receptor_PDB")
  if (!dir.exists(ligand_dir)) dir.create(ligand_dir, recursive = TRUE)
  if (!dir.exists(receptor_dir)) dir.create(receptor_dir, recursive = TRUE)

  # ---- Download protein structure ----
  # always download ligands (since docking_ligand_dir may depend on them)
  download_structures(ligand_match, ligand_dir, "ligand")

  # only download receptors if user did NOT provide docking_receptor_dir
  if (is.null(docking_receptor_dir)) {
    download_structures(receptor_match, receptor_dir, "receptor")
  }

  # ---- Ligand structure source ----
  cas_dir <- NULL

  if (use_fda) {
    # ---- FDA compounds support ----
    if (is.null(fda_txt_path) || !file.exists(fda_txt_path)) {
      stop("[FDA] use_fda = TRUE but fda_txt_path is missing or does not exist!")
    }

    fda_dir <- file.path(output_dir, "fda_results")
    if (!dir.exists(fda_dir)) dir.create(fda_dir, recursive = TRUE)

    fda_cas <- trimws(readLines(fda_txt_path))
    if (length(fda_cas) > 0) {
      message("[FDA] Downloading structures for ", length(fda_cas), " FDA compounds...")
      system2(python_bin, args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        fda_txt_path, fda_dir
      ))
      docking_ligand_dir <- fda_dir
    } else {
      stop("[FDA] fda.txt is empty! Please provide CAS numbers.")
    }

  } else {
    # ---- CAS compounds support ----
    if (!is.null(cas_txt_file) && file.exists(cas_txt_file)) {
      cas_dir <- file.path(output_dir, "ligand_structure")
      if (!dir.exists(cas_dir)) dir.create(cas_dir)
      cas_numbers <- trimws(readLines(cas_txt_file))
      message("[CAS] Downloading structures for ", length(cas_numbers), " compounds...")
      system2(python_bin, args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        cas_txt_file, cas_dir
      ))
      docking_ligand_dir <- cas_dir
    }
  }

  # ---- AutoDock Vina docking ----
  docking_ligand_dir <- if (is.null(docking_ligand_dir)) cas_dir else docking_ligand_dir

  if (is.null(docking_receptor_dir)) {
    docking_receptor_dir <- list(ligand_dir, receptor_dir)
  } else {
    message("[INFO] Using user-provided receptor directory: ", docking_receptor_dir)
  }

  receptor_files <- unlist(lapply(docking_receptor_dir, function(d) {
    list.files(d, pattern = "_prepared\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  }))

  ligand_files <- list.files(docking_ligand_dir, pattern = "\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  ligand_files <- ligand_files[file.size(ligand_files) > 0]

  message("Starting docking...")
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
  }
  message("✅ Docking complete.")
}
