# Use as follows:
# python aws_upload.py <full path to your AWS creds sh file> <full path to your local file to upload> <S3 URI to upload to, e.g. s3://immuneaging/scratch>

import os
import sys
from utils import set_access_keys

if __name__ == "__main__":
    creds_file = sys.argv[1]
    fn = sys.argv[2].split("/")[-1]
    dn = sys.argv[2][:-len(fn)]
    s3_path = sys.argv[3]
    set_access_keys(creds_file)
    sync_cmd = 'aws s3 sync {} {} --exclude "*" --include {}'.format(dn, s3_path, fn)
    os.system(sync_cmd)
