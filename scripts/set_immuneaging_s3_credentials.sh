#! /bin.bash

#####################
# DO NOT DISTRIBUTE #
#####################

# Sets AWS credentials for First Name, Last Name

# usage: bash set_immuneaging_s3_credentials.sh

# This script will set your personal access keys for the Yosef Lab's immuneaging s3 bucket.
# This script should be run once per bash session and will only allow access in said
# bash session. Once the session is closed, access will cease and the aws credentials will 
# revert back to the aws config credentials.

export AWS_ACCESS_KEY_ID=<KEY>
export AWS_SECRET_ACCESS_KEY=<KEY>
