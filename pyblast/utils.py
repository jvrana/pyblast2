import shlex
import subprocess


def run_cmd_str(cmd_str):
    print("CMD: "+cmd_str)
    args = shlex.split(cmd_str)
    output = subprocess.Popen(args)
    output.wait()
    return output


def run_cmd(cmd, **kwargs):
    run_cmd_str(dict_to_cmd(cmd, **kwargs))


def dict_to_cmd(cmd, **kwargs):
    return cmd+" "+' '.join(["-{} {}".format(k, kwargs[k]) for k in kwargs])


def str_to_f_to_i(v):
    try:
        v = float(v)
    except ValueError as e:
        pass
    try:
        v = int(v)
    except ValueError as e:
        pass
    return v
