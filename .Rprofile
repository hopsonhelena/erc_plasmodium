if (!nzchar(Sys.getenv("RENV_PATHS_ROOT"))) {
  Sys.setenv(RENV_PATHS_ROOT = file.path(getwd(), "renv", ".cache"))
}
if (!nzchar(Sys.getenv("RENV_PATHS_CACHE"))) {
  Sys.setenv(RENV_PATHS_CACHE = file.path(getwd(), "renv", ".cache", "cache"))
}
# CHPC's shared R installation does not need renv's copied system-library sandbox.
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
source("renv/activate.R")
