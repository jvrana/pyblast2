"""cmd.py"""

import shlex
from subprocess import check_output, STDOUT, CalledProcessError


def run_cmd_str(cmd_str):
    """Runs a command from a string"""
    print("CMD: " + cmd_str)
    args = shlex.split(cmd_str)
    try:
        output = check_output(args, stderr=STDOUT)
    except CalledProcessError as error:
        error_message = ">>> Error while executing:\n{}\n>>> Returned with error:\n{}".format(cmd_str, str(error.output))
        error_message += "\nParameter info: https://www.ncbi.nlm.nih.gov/books/NBK279684/"
        raise Exception(error_message)

    # output = subprocess.Popen(args, shell=False)
    # output.wait()
    return output.decode()


def run_cmd(cmd, **kwargs):
    """Run a command using parameters kwargs"""
    run_cmd_str(dict_to_cmd(cmd, **kwargs))


def dict_to_cmd(cmd, **kwargs):
    """Create a command string for cmd and parameters 'kwargs'"""
    return cmd + " " + ' '.join(["-{} {}".format(k, kwargs[k]) for k in kwargs])

