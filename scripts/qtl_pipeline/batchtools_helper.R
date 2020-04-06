#' Extension to the waitForJobs function of the BatchJobs package which shows
#' some strange behaviour when waiting for jobs (database locked)
#' so we need to make it extra failsafe.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param reg batchtools registry
#' @param waittime time to wait before updating job status
#' @param nretry number of time to retry getting the job status before throwing
#'        an error
#'
#' @imports batchtools
#' 
#' 

# rfun = do.stuff
# idx = list(1:3)
# more.args = list()
# name = "test"
# dir = "tmp_dir_test"
# clean.up = F
# resources = list(time="00:05:00", memory="200M")

myWaitForJobs <- function(reg, waittime=3, nretry=100) {
  require(batchtools)
  success = FALSE
  while (nretry > 0 && !success) {
    status = tryCatch({
      while (TRUE) {
        status = getStatus(reg = reg)
        if (status$done + status$error + status$expired == status$defined) {
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


run.batchtools <- function(rfun, idx, more.args, name, dir, clean.up=TRUE,
                           resources=list(), n.chunks=NULL, reuse.registry=FALSE) {
  require(batchtools)
  
  if (reuse.registry) {
    reg <- readRDS(paste0(dir, "/registry.rds"))
  } else {
    reg = makeRegistry(file.dir=dir)
  }
  print("Registry created")
  batchMap(fun = rfun,
           idx,
           more.args=more.args,
           reg = reg)
  
  # ## jobs can be packed together in chunks if there are too many
  # chunked = findJobs(reg = reg)
  # if (!is.null(n.chunks)) {
  #   chunked = chunk(chunked, n.chunks=n.chunks, shuffle = TRUE)
  # }
  # print("checked for chunks")
  
  ## job delay to prevent concurrent access to the database by too many jobs
  print("Let's submit jobs")
  submitJobs(reg = reg, resources=resources)
  print("jobs submitted")
  Sys.sleep(20)
  
  print("Let's wait!")
  
  ## also a custom wait function that is more error tolerant with many jobs
  myWaitForJobs(reg, waittime=30, nretry=100)
  Sys.sleep(20)
  res = reduceResultsList(reg=reg)
  print(head(getJobTable(reg=reg), 10))
  if (clean.up) {
    removeRegistry(wait=0, reg=reg)
  }
  return(res)
}