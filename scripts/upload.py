import pandas as pd
import argparse
import warnings
import os


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


def read_immune_aging_sheet(sheet, output_fn=None, sheet_name=None):
    url = "https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/gviz/tq?tqx=out:csv&sheet={}".format(
        sheet
    )

    # testing url of bad sample sheet
    # url = "https://docs.google.com/spreadsheets/d/1YO1HLGLnO3PPUiK1vKZd52yoCInwpLl60zoi4zxkOrE/gviz/tq?tqx=out:csv&sheet={}".format(
    #     sheet
    # )

    try:
        import gdown
    except ImportError as e:
        raise ImportError(
            "gdown is not installed. Please install gdown via: pip install gdown"
        )

    with warnings.catch_warnings(record=True) as w:
        warnings.filterwarnings("error")
        output_fn = gdown.download(url, output_fn, quiet=False)

        if len(w) == 1:
            warnings.showwarning(
                msg.message, msg.category, msg.filename, msg.lineno, msg.line
            )

    data = pd.read_csv(output_fn)
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


def get_fastq_gzs_in_folder(folder, recursive=False):
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
        files = os.listdir(folder)
    files = [f for f in files if f.endswith(".fastq.gz")]
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
        "--destination", help="destination folder to upload data to", required=True
    )
    parser.add_argument(
        "--fastq",
        help="path to folder containing fastq.gz files",
        required=True,
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        default=False,
        help="if flag set, won't check subfolders",
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


def check_sheet(donors, samples, dictionary):
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
            v = str(v)
            # for numerical ranges
            if "-" in v:
                tmp_v = v.split("-")
                for x in tmp_v:
                    try:
                        float(x)
                    except:
                        msg = "The range {} is not numeric. Please edit the online Google sheet and rerun script.".format(
                            tmp_v
                        )
                        warnings.warn(msg)
                        is_valid = False
            else:
                try:
                    float(v)
                except:
                    msg = "{} in {} is not numeric. Please edit the online Google sheet and rerun script.".format(
                        v, col_name
                    )
                    warnings.warn(msg)
                    is_valid = False
        return is_valid

    is_valid = True

    dictionary["GEX chem"] = dictionary["Chemistry"]
    dictionary["CITE chem"] = dictionary["Chemistry"]
    dictionary["TCR chem"] = dictionary["Chemistry"]
    dictionary["BCR chem"] = dictionary["Chemistry"]
    del dictionary["Chemistry"]

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
    height_is_valid = check_col_is_numerical(donors, "height (cm)")
    smoker_is_valid = check_col_is_numerical(donors, "smoker (pack-years)")

    is_valid = (
        is_valid
        and bmi_is_valid
        and age_is_valid
        and height_is_valid
        and smoker_is_valid
    )

    if not is_valid:
        raise ValueError(
            "Invalid values in metadata sheet. Check above warnings. "
            #"Please fix and rerun command or use -f to force upload."
            "Please fix and rerun command."
        )

    return is_valid


def check_fastq_filenames(fastq_fns, samples, donors):
    """
    Check that the fastq filenames conforms to data schema standards:
    <donor id>_<library_type>_<library ID>_<10x reaction>_<seq run>_<lane>_<read>.fastq.gz

    Parameters
    ----------
    fastq_fns
        filenames of all fastqs to be uploaded
    """
    fn_format = "<donor_id>_<library_type>_<library_ID>_<10x_reaction>_<seq_run>_<lane>_<read>.fastq.gz"
    
    fns = [f.split("/")[-1] for f in fastq_fns]
    donors.index = donors["Donor ID"]
    donor_ids = list(donors["Donor ID"])

    samples.index = samples["Donor ID"]

    # ADT is a CITE-seq library, often referred to as ADT or antibody-derived tag.
    lib_types = ["GEX", "ADT", "TCR", "BCR"]
    is_valid = True

    for f in fns:
        try:
            data = f.split(".")
            assert ".".join(data[-2:]) == "fastq.gz"
            data = data[0].split("_")
            donor_id = data[0]
            lib_type = data[1]
            lib_id = data[2]
            tenx_reaction = data[3]
            seq_run = data[4]
            lane = data[5]
            read_num = data[6]
        except ImportError as e:
            raise ValueError("Bad file name. Please adhere to the format " + fn_format)

        # check donor id
        valid_donor_id = donor_id in donor_ids

        # check library type
        valid_lib_type = lib_type in lib_types

        # check library id
        valid_lib_id = True
        if valid_donor_id and valid_lib_type:

            if samples.loc[donor_id][lib_type + " lib"].isnull().values.any():
                msg = "{} has no entry in the '{} lib' column in the Samples".format(
                    f, lib_type
                )
                warnings.warn(msg)
                valid_lib_id = False
            else:
                l_ids = samples.loc[donor_id][lib_type + " lib"]

                # sample id is not unique so get all library id values
                if isinstance(l_ids, pd.core.series.Series):
                    l_ids = l_ids.values
                elif isinstance(l_ids, str):
                    l_ids = [l_ids]

                # some library ids are seperated by commas
                l_ids = [s.split(",") for s in l_ids]

                # if they were seperated by commas, collapse the list
                if isinstance(l_ids[0], list):
                    l_ids = [item for sublist in l_ids for item in sublist]

                l_ids = [
                    s.strip() for s in l_ids
                ]  # l_ids here is finally all library ids for the samples of a donor

                valid_lib_id = lib_id in l_ids

        else:
            valid_lib_id = False


        valid_tenx_reaction = False
        valid_seq_run = False
        valid_lane = False
        valid_read_num = False
        
        try:
            int(valid_tenx_reaction)
            valid_tenx_reaction = True
        except: pass

        try:
            int(seq_run)
            valid_seq_run = True
        except: pass
            
        try:
            assert lane[0] == "L"
            int(lane[1:])
            valid_lane = True
        except: pass
            
        try:
            assert read_num[0] == "R"
            int(read_num[1:])
            valid_read_num = True
        except: pass

        valid_sample = valid_donor_id and valid_lib_type and valid_lib_id and valid_tenx_reaction and valid_seq_run and valid_lane and valid_read_num

        if valid_sample is False:
            is_valid = False
            msg = "Error for file: {}. ".format(f)
            if valid_donor_id is False:
                msg += "Donor id does not exist in Sample Sheet. "
            if valid_lib_type is False:
                msg += "Library type of '{}' is invalid. Options: ['GEX', 'ADT', 'TCR', 'BCR'] .".format(
                    lib_type
                )
            if valid_donor_id and valid_lib_type and (valid_lib_id is False):
                msg += "Library id of {} is not in Sample Sheet.".format(
                    lib_id)
            if valid_tenx_reaction is False:
                msg += "10x_reaction number in the file name should be a number (e.g., 0001); instead found {}. ".format(tenx_reaction)
            if valid_seq_run is False:
                msg += "seq_run number in the file name should be a number (e.g., 001); instead found {}. ".format(seq_run)
            if valid_lane is False:
                msg += "lane in the file name should be in a valid format (e.g., L001); instead found {}. ".format(lane)
            if valid_read_num is False:
                msg += "read in the file name should be in a valid format (e.g., R1); instead found {}. ".format(read_num)
            if (valid_tenx_reaction and valid_seq_run and valid_lane and valid_read_num) is False:
                msg += "Please adhere to the format {}. ".format(fn_format)
            warnings.warn(msg)

    if is_valid:
        print("fastq filename check successful.")
    else:
        raise ValueError("Error in fastq filenames. Check above warnings.")

    return is_valid


def upload_to_s3(source, destination):
    """
    Upload source filenames to destination in s3.
    """
    if destination == "sanger" or destination == "columbia":
        destination = "raw_" + destination
    else:
        destination = "test_folder"
    aws_cmd = "aws s3 sync {} s3://immuneaging/{}".format(source, destination)
    print("Running aws command: ", aws_cmd)
    os.system(aws_cmd)


def validate_args(args):
    destinations = ["test", "sanger", "columbia"]
    print(args)
    if args.destination not in destinations:
        raise ValueError(
            "--destination must be one of the followings: {0}".format(", ".join(destinations)))


if __name__ == "__main__":
    args = parse_args()
    validate_args(args)

    fastq_fns = get_fastq_gzs_in_folder(args.fastq, recursive=args.recursive)

    # comment out if we want to give the force option
    args.force = False

    if not args.force:
        dict_sheet = read_immune_aging_sheet(sheet="Dictionary")
        dictionary = make_immuneaging_dictionary(dict_sheet)
        donors = read_immune_aging_sheet(sheet="Donors")
        samples = read_immune_aging_sheet(sheet="Samples")
        check_sheet(donors, samples, dictionary)

        check_fastq_filenames(fastq_fns, samples, donors)

    set_access_keys(args.aws_keys)
    upload_to_s3(args.fastq, args.destination)
