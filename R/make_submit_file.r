#' @title Submit batch jobs based on batch specification
#'
#' @description This function writes a bash script that submits
#'              batch jobs to the SLURM scheduler based
#'              on the output of `make_batch_script()`.
#'
#' @param specname A character string that should be the same as the
#'                 `specname` argument provided to `make_batch_script()`.
#' @param specdir A character string giving the full path to the destination
#'                directory for the script.
#' @param batches A list containing vectors of integers, where the position
#'                in the list corresponds to the batch ID.
#' @param offsets A numeric vector giving the offsets for each batch.
#'
#' @export
make_submit_file <- function(specname, specdir, batches, offsets) {
  submit_file <- file.path(specdir, paste0(specname, "_batch_submit.sh"))
  if (file.exists(submit_file)) {
    file.remove(submit_file)
  }

  # create the submit file
  for (i in seq_along(batches)) {
    sink(submit_file, append = TRUE)
    cat("sed 's/$OFFSET/", offsets[i], "/' ",
        paste0(specname, "_batch_template.sh | sbatch"),
        "\n", sep = "")
    sink()
  }
}
