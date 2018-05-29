"""cmd.py"""

import shlex
from subprocess import check_output


def run_cmd_str(cmd_str):
    """Runs a command from a string"""
    print("CMD: " + cmd_str)
    args = shlex.split(cmd_str)
    output = check_output(args)
    # output = subprocess.Popen(args, shell=False)
    # output.wait()
    return output.decode()


def run_cmd(cmd, **kwargs):
    """Run a command using parameters kwargs"""
    run_cmd_str(dict_to_cmd(cmd, **kwargs))


def dict_to_cmd(cmd, **kwargs):
    """Create a command string for cmd and parameters 'kwargs'"""
    return cmd + " " + ' '.join(["-{} {}".format(k, kwargs[k]) for k in kwargs])

