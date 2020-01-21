# generate_ms_call

# Simple ms call ----------------------------------------------------------

gc <- 5:10 # gene copies or "n"
mu <- 1e-8 # mutation rate
Ne <- 10 # effect pop size (could be a vector)
r <- 2*mu # scaled rate

popsize <- 10 # true popsize (this is a vector)

theta <- 4*Ne*mu*popsize # theta

options(scipen = 999)

# now generate the actual call
call <- paste0("ms ", gc, " ", popsize, " -t ", theta)
call # double check

# Generate PARAMS ONLY
# now generate the actual call
#call <- paste0(gc, " ", popsize, " -t ", theta)
#call # double check


# write out
readr::write_lines(call, path = "data/ms_params_tst_001.txt")

# append to existing file?
readr::write_lines(call, path = "data/ms_params_tst_001.txt", append = T)

#cat(call, file="output.txt",sep="\n",append=TRUE)

