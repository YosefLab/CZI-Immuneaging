# conda activate immune_aging.py_env.v2
# Run: Rscript launch_vision.R /path/to/vision_object.rds port_number
# the session will be accessible via http://<node>.cs.berkeley.edu:<port>/Results.html where <node> is the node in which the session was executed and <port> is the requested port (assuming the port is open).
# e.g., port 8470 is open on s124 or s125
library("VISION")
args = commandArgs(trailingOnly=TRUE)
input_file <- args[1]
port <- args[2]
session_name <- args[3]
vis <- readRDS(input_file)
viewResults(vis, host='0.0.0.0', browser=FALSE, port = as.numeric(port), name = session_name)