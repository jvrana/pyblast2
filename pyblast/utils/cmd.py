"""cmd.py."""
import shlex
from subprocess import CalledProcessError
from subprocess import check_output
from subprocess import PIPE
from subprocess import Popen
from subprocess import STDOUT


def run_cmd_str(cmd_str):
    """Runs a command from a string."""
    args = shlex.split(cmd_str)
    try:
        process = Popen(args, shell=False, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as error:
        error_message = ">>> Error while executing:\n{}\n>>> Returned with error:\n{}".format(
            cmd_str, str(error.output)
        )
        error_message += (
            "\nParameter info: https://www.ncbi.nlm.nih.gov/books/NBK279684/"
        )
        raise Exception(error_message)

    (stdoutdata, stderrdata) = process.communicate(timeout=120)
    if stderrdata:
        raise CalledProcessError(stderrdata)
    return stdoutdata.decode()


def run_cmd(cmd, **kwargs):
    """Run a command using parameters kwargs."""

    flags = [k for k, v in kwargs.items() if v is None]
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    run_cmd_str(dict_to_cmd(cmd, flags, **kwargs))


def dict_to_cmd(cmd, flags, **kwargs):
    """Create a command string for cmd and parameters 'kwargs'."""
    cmd_str = cmd + " " + " ".join(["-{} {}".format(k, kwargs[k]) for k in kwargs])
    if flags:
        flag_str = " ".join(["-" + f for f in flags])
    else:
        flag_str = ""
    return cmd_str + " " + flag_str
