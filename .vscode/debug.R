env <- Sys.getenv()
envnames <- names(env)
rnames <- envnames[startsWith(envnames, "R_")]
cached_names <- rnames
ld_lib_path <- Sys.getenv("LD_LIBRARY_PATH")
if (ld_lib_path != "") {
    cached_names <- c("LD_LIBRARY_PATH", rnames)
}
writeLines(paste0(cached_names, "=", env[cached_names]), ".vscode/.env")