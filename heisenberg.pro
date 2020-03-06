; command line calling file for calls to HEISENBERG using iterative Fourier diffuse filtering

if keyword_set(COMMAND_LINE_ARGS()) then master_inputfile=COMMAND_LINE_ARGS() $
                                    else read,' please specify the full/absolute path of the input file (enter 0 to stop autorun): ',master_inputfile
.com diffuse_iteration

diffuse_iteration, master_inputfile
