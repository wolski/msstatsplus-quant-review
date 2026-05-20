## timing.R - two-phase wall-clock timing shared by all models_*.R.

tic <- function() proc.time()[3]
toc <- function(t0) as.numeric(proc.time()[3] - t0)

write_timing <- function(method, out_path, preprocess_seconds, model_seconds) {
  data.table::fwrite(
    data.frame(method = method,
               preprocess_seconds = preprocess_seconds,
               model_seconds      = model_seconds),
    file = file.path(out_path, paste0(method, "_timing.csv"))
  )
  invisible(NULL)
}
