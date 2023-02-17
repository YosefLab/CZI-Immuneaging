import os

def set_access_keys(filepath, return_dict = False):
    """
    Sets the user's access keys to the AWS S3 bucket.

    Assumes this file includes only two uncommented lines with keys:
    export AWS_ACCESS_KEY_ID=<key>
    export AWS_SECRET_ACCESS_KEY=<key>

    If return_dict == True then only returns the dictionary.
    """
    keys = {}
    with open(filepath) as fp:
        for i in fp:
            if len(i.rstrip()) and i[0]!="#":
                cmd = i.rstrip().split(" ")
                assert(cmd[0] == "export")
                cmd.pop(0)
                for i in range(len(cmd)-1,-1,-1):
                    if len(cmd[i]) == 0:
                        cmd.pop(i)
                pair = "".join(cmd).split("=")
                keys[pair[0]] = pair[1]
    if return_dict:
        return keys
    for k in keys:
        os.environ[k] = keys[k]