# load cluster functions
cluster.functions <- makeClusterFunctionsSlurm(template="~/ext_tools/batchtools/template_slurm.tmpl",
                                               array.jobs=TRUE)
mail.start = "none"
mail.done = "none"
mail.error = "none"

# set some default resources
default.resources <- list(partition="my_queue",
                          memory="8G",
                          ncpus = 1,
                          measure.memory = TRUE,
                          walltime="12:00:00")