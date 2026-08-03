# Source all project functions from any script in scripts/.
siga_project_root <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command, value = TRUE)
  if (length(file_argument)) {
    script <- normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = FALSE)
    return(dirname(dirname(script)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SIGA_ROOT <- siga_project_root()
source(file.path(SIGA_ROOT, "R", "siga_core.R"), local = FALSE)
source(file.path(SIGA_ROOT, "R", "siga_study_definitions.R"), local = FALSE)
source(file.path(SIGA_ROOT, "R", "siga_workflows.R"), local = FALSE)
source(file.path(SIGA_ROOT, "config", "manuscript_settings.R"), local = FALSE)
