#' Extension to the waitForJobs function of the BatchJobs package which shows
#' some strange behaviour when waiting for jobs (database locked)
#' so we need to make it extra failsafe.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param reg BatchJobs registry
#' @param waittime time to wait before updating job status
#' @param nretry number of time to retry getting the job status before throwing
#'        an error
#'
#' @imports BatchJobs
myWaitForJobs <- function(reg, waittime=3, nretry=100) {
  require(BatchJobs)
  success = FALSE
  while (nretry > 0 && !success) {
    status = tryCatch({
      while (TRUE) {
        status = showStatus(reg)
        if (status$done + status$expired == status$n) {
          cat("done\n")
          return(list(success=TRUE, nretry=nretry))
        }
        Sys.sleep(waittime)
      }
      return(list(success=FALSE, nretry=nretry))
    }, error=function(e) {
      cat("Error while waiting for jobs:\n")
      print(e)
      cat("\nnumber of retries left: ", nretry - 1, "\n")
      Sys.sleep(waittime + runif(1, 0, 3))
      return(list(success=FALSE, nretry=nretry - 1))
    })
    success = status$success
    nretry = status$nretry
    cat("success after the tryCatch block:", success, "\n")
    cat("nretry after the tryCatch block:", nretry, "\n")
  }
  
  
  if (!success) {
    err.msg = paste("Error during batch processing in registry")
    save(envir=sys.frame(), list=ls(envir=sys.frame()), file=file.path(dir, "error_image.RData"))
    stop(err.msg)
  }
}


run.batchjobs <- function(rfun, idx, more.args, name, dir, clean.up=TRUE, resources=list(), n.chunks=NULL, reuse.registry=FALSE) {
  require(BatchJobs)
  
  if (reuse.registry) {
    load(paste0(dir, "/registry.RData"))
  } else {
    reg = makeRegistry(name, file.dir=dir)
  }
  batchMap(reg, rfun, idx, more.args=more.args)
  
  ## jobs can be packed together in chunks if there are too many
  chunked = getJobIds(reg)
  if (!is.null(n.chunks)) {
    chunked = chunk(chunked, n.chunks=n.chunks, shuffle = TRUE)
  }
  
  ## job delay to prevent concurrent access to the database by too many jobs
  submitJobs(reg, ids=chunked, resources=resources, job.delay=function(n, i) runif(1, 0, 0.1 * n))
  Sys.sleep(20)
  
  ## also a custom wait function that is more error tolerant with many jobs
  myWaitForJobs(reg, waittime=30, nretry=100)
  res = reduceResultsList(reg)
  if (clean.up) {
    removeRegistry(reg, ask="no")
  }
  return(res)
}


## Mini example:
#
#
# rfun <- function(i, filename) {
#   # do some very important stuff
#   return(i)
# }
# 
# test_batch <- run.batchjobs(rfun, 1:10, more.args=list(filename=123),
#                             "simulation", "tmp_dir_simulation")
# saveRDS(test, "test_batch_jobs.RDS")
#
# load("/home/icb/ines.assum/rstudio/tmp_dir_simulation/registry.RData")
# library(BatchJobs, lib.loc = "/mnt/storageGluster/groups/groups_epigenereg/packages/2017/R/3.4")
# submitJobs(reg, findErrors(reg))
# submitJobs(reg, findExpired(reg))