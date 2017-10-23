import shlex
import subprocess
import os

# TODO: bin/sh is broken in Travis-CI linux >>>  E OSError: [Errno 8] Exec format error: 'makeblastdb'

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def run_cmd_str(cmd_str):
    print("CMD: "+cmd_str)
    args = shlex.split(cmd_str)
    output = subprocess.Popen(args, shell=False)
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
