import pandas as pd
import argparse
import warnings


def set_access_keys(filepath):
    """
    Sets the user's access keys to the AWS S3 bucket.

    Assumes this file includes only two uncommented lines with keys:
    export AWS_ACCESS_KEY_ID=<key>
    export AWS_SECRET_ACCESS_KEY=<key>
    """
    keys = {}
    with open(filepath) as fp:
        for i in fp:
            if len(i.rstrip()) and i[0] != "#":
                cmd = i.rstrip().split(" ")
                assert cmd[0] == "export"
                cmd.pop(0)
                for i in range(len(cmd) - 1, -1, -1):
                    if len(cmd[i]) == 0:
                        cmd.pop(i)
                pair = "".join(cmd).split("=")
                keys[pair[0]] = pair[1]
    for k in keys:
        os.environ[k] = keys[k]


def read_google_sheet(url, output_fn=None, sheet_name=None):
    """
    Downloads and reads a Google Sheet. Must be xlsx.

    Parameters
    ----------
    url
        url of publicly accessible Google Sheet.
    output_fn
        Optional output filename
    sheet_name
        Specific sheet to return. By default returns all sheets in spreadsheet.
    """
    try:
        import gdown
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "gdown is not installed. Please install gdown via: pip install gdown"
        )

    # url needs to be in a specific format
    with warnings.catch_warnings(record=True) as w:
        warnings.filterwarnings("error")
        output_fn = gdown.download(url, output_fn, quiet=False)
        if len(w) == 1:
            warnings.showwarning(
                msg.message, msg.category, msg.filename, msg.lineno, msg.line
            )

    # read spreadsheet
    data = pd.read_excel(output_fn, sheet_name=None)
    return data


def make_immuneaging_dictionary(df):
    """
    Makes the dictionary from immuneaging dictionary

    Parameters
    ----------
    df
        "Dictionary" sheet from Metadata Google Sheet as pd.DataFrame
    """
    # Ignore empty columns
    col_names = [c for c in df.columns if not c.startswith("Unnamed:")]

    # following line gets the non_null values for each column of the dataframe
    # then makes it into a dictionary where key is the column name and value is list of nonnull values
    dictionary = {k: list(df[df[k].notnull()][k]) for k in col_names}
    return dictionary


def get_fastq_gz_in_folder(folder, recursive=False):
    """
    Get all the fastqs in a folder

    Parameters
    ----------
    folder
        folder to get the fastq.gz from
    recursive
        if True, will check subfolders as well
    """
    if recursive:
        import glob

        files = glob.glob(folder + "/**/*.fastq.gz", recursive=True)
    else:
        import os

        all_dir = os.listdir(folder)
    files = [f for f in all_dir if f.endswith(".fastq.gz")]
    return files


def parse_args():
    """
    Parse command line arguments

    Running `python upload.py -h` will show arguments and how to use.
    """
    parser = argparse.ArgumentParser(
        description="Upload files to immuneaging s3 bucket."
    )
    parser.add_argument("--aws_keys", help="path to aws keys", required=True)
    parser.add_argument(
        "--fastq",
        help="path to folder containing fastqs",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        default=False,
        help="if flag set, will check all subfolders",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="if force upload, ignore checks",
    )
    args = parser.parse_args()
    return args


def get_non_null_values(df, col_name):
    return df[df[col_name].notnull()][col_name]


def check_sheet(metadata):
    """
    Checks Metadata Sheet.

    Parameters
    ----------
    url
        url of publicly accessible Metadata Google Sheet.
    """
    def check_col_is_numerical(df, col_name):
        vals = get_non_null_values(df, col_name)
        is_valid = True
        for v in vals:
            if not str(v).isnumeric():
                msg = "{} in {} is not numeric. Please edit the online Google sheet and rerun script.".format(
                    v, col_name
                )
                warnings.warn(msg)
                is_valid = False
        return is_valid

    dictionary = make_immuneaging_dictionary(metadata_df["Dictionary"])

    donors = metadata_df["Donors"]
    samples = metadata_df["Samples"]

    is_valid = True
    for col_name, valid_values in dictionary.items():
        if col_name in donors.keys():
            values = get_non_null_values(donors, col_name)
        elif col_name in samples.keys():
            values = get_non_null_values(samples, col_name)
        else:
            msg = "{} not in Donors or Samples page.".format(col_name)
            warnings.warn(msg)
            values = []

        for val in values:
            if val in valid_values:
                pass
            else:
                warnings.warn(
                    "{} not a valid value in column {}. "
                    "Please edit the online Google sheet and rerun script".format(
                        val, col_name
                    )
                )
                is_valid = False

    bmi_is_valid = check_col_is_numerical(donors, "BMI (kg/m^2)")
    age_is_valid = check_col_is_numerical(donors, "Age (years)")

    is_valid = (is_valid and bmi_is_valid and age_is_valid)

    if not is_valid:
        raise ValueError(
            "Invalid values in metadata sheet. Check above warnings. "
            "Please fix and rerun command or use -f to force upload."
        )

    return metadata_df

def check_fastq_filenames(fastq_fns, metadata_df):
    """
    Check that the fastq filenames conforms to data schema standards

    Parameters
    ----------
    fastq_fns
        filenames of all fastqs to be uploaded
    metadata_df
        DataFrame containing metadata info

    """
    raise NotImplementedError

def upload_to_s3(source_fns, destination):
    """
    Upload source filenames to destination in s3.
    """
    raise NotImplementedError

if __name__ == "__main__":
    args = parse_args()

    # set_access_keys(args.aws_keys)
    # fastq_fns = get_fastq_gz_in_folder(args.fastq, recursive=args.recursive)

    url = "https://drive.google.com/uc?id=1dzd6WjPaki1plsG2phgcAnvXjfbrNDkR"
    metadata_df = read_google_sheet(url)

    if not args.force:
        check_sheet(metadata_df)
        # check_fastq_filenames(fastq_fns, metadata_df)

    set_access_keys(args.aws_keys)
    upload_to_s3()
