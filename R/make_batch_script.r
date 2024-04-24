#' @title Make a batch script to submit a SLURM job
#' @description This function writes a batch script that contains
#'              specifications for a job array that will run the
#'              requested number of simulations.
#'
#' @param specname A character string giving the name of the job specification.
#'                 This parameter will be used to direct logs and output to the
#'                 correct locations.
#' @param walltime A character string giving the walltime for the job.
#' @param partition A character string giving the partition for the job.
#' @param mem A character string giving the memory for the job.
#' @param ncores An integer giving the number of cores for the job.
#' @param array A character string giving the array ID range.
#' @param offset A character string giving the offset that will map the array
#'               IDs to the correct simulation number. Environment variable.
#' @param logname A character string giving the naming format of the log file.
#'                Log files will be written to the .logs directory within the
#'                directory specified by specname.
#' @param batchid An integer giving the batch id for the job.
#' @param nsims An integer giving the number of simulations runs for the job.
#'              Usually, this is set to 1 because we run 1 simulation per array
#'              ID. Environment variable.
#' @param nsteps An integer giving the number of time steps in the simulation.
#'               Environment variable.
#' @param add_arrivals A numeric input indicating the number of arrivals to add
#'                     per timestemp. Chosen in preliminary simulations to
#'                     maintain the population size at approximately 20,000.
#'                     Environment variable.
#' @param rootdir Full path to the project root directory.
#' @param simdir Full path to the directory in which to write the
#'               simulation output. This should typically be a directory within
#'               the ~/scratch directory.
#' @param rscript_exec Full path to the `Rscript` executable.
#' @param rscript_exec_args Arguments to pass to the `Rscript` executable.
#'                          Defaults to "--vanilla".
#' @param rscript_file Path to the R script to execute relative to project `rootdir`.
#'                     Will be passed as a full path to the SLURM scheduler.
#' @param make_submit Boolean. If `TRUE`, will run `make_submit_file()`.
#' @param ... Pass arguments to `make_submit_file()`.
#'
#' @return A batch script that contains specifications for a SLURM job array.
#'         This specification is referred to as a template and is used as the
#'         basis for the job submission file written by `make_submit_file()`.
#'
#' @importFrom here here
#' @export
make_batch_script_template <- function(specname, walltime, partition, mem,
                                       ncores, array, offset, logname,
                                       nsims, nsteps, add_arrivals,
                                       rootdir, simdir,
                                       rscript_exec = "/users/jgantenb/build-sources/R-4.3.2/bin/Rscript",
                                       rscript_exec_args = "--vanilla",
                                       rscript_file,
                                       make_submit = TRUE,
                                       ...) {

  # set up specification directory (create if doesn't exist)
  specdir <- here("simulator_run", specname)
  if (!dir.exists(specdir)) {
    dir.create(specdir)
  }

  # set up specification .logs directory (create if it doesn't exist)
  if (!dir.exists(file.path(specdir, ".logs"))) {
    dir.create(file.path(specdir, ".logs"))
  }

  # full path to the log file for each simulation
  log_fullpath <- file.path(specdir, ".logs", logname)

  sb <- "#SBATCH"

  specs <- paste(
    "#!/bin/bash",
    paste(sb, "-J", specname),
    paste0(sb, " --time=", walltime),
    paste(sb, "-p", partition),
    paste0(sb, " --mem=", mem),
    paste(sb, "-n", ncores),
    paste0(sb, " --array=", array),
    paste(sb, "-o", log_fullpath),
    paste0(
      sb,
      " --export=ALL,NSIMS=", nsims,
      ",NSTEPS=", nsteps,
      ",ARRIVE_RATE_ADD_PER20K=", add_arrivals,
      ",SIMDIR=", simdir,
      ",SLURM_ARRAY_MAP_OFFSET=", offset
    ),
    paste(sb, "--mail-type=ALL"   ),
    paste(sb, "--mail-user=jrgant@brown.edu"),
    paste("cd", rootdir, "|| return"),
    paste(rscript_exec, paste0(rootdir, "/", rscript_file), rscript_exec_args),
    sep = "\n"
  )

  writeLines(
    text = specs,
    con = file.path(
      specdir,
      paste0(specname, "_batch_template.sh")
    )
  )

  if (make_submit) {
    make_submit_file(
      specname = specname,
      specdir = specdir,
      ...
    )
  }

}
