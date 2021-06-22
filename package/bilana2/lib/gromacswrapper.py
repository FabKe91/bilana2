'''
    API to easily run gromacs commands

'''
import os
import sys
import subprocess
from datetime import datetime

def find_executable(executable, path=None):
    """
    FROM:
    # https://gist.github.com/4368898
    # Public domain code by anatoly techtonik <techtonik@gmail.com>
    # AKA Linux `which` and Windows `where`

    Find if 'executable' can be run. Looks for it in 'path'
    (string that lists directories separated by 'os.pathsep';
    defaults to os.environ['PATH']). Checks for all executable
    extensions. Returns full path or None if no command is found.
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (_, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (_, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    raise LookupError("No Gromacs executable found")

def exec_gromacs(gmx, gmxlib, inpargs, interactive=None, backup=True):
    '''
        Execute Gromacs commands.
        cmd variable stores the actual terminal command (can be anything,
        but here is specifically for gromacs)
        cmd must be a list with items being the strings (no whitespaces) of the command
            "['gmx cmd', '-f', 'tprfile', '-e', 'en.edr']"
        if a command needs user input (STDIN, like in gmx make_ndx) inp_str can be parsed
        as string
    '''
    if backup:
        backup_flag = "-backup"
    else:
        backup_flag = "-nobackup"

    if isinstance(inpargs, list):
        cmd = [gmx, backup_flag, "-nocopyright", gmxlib, *inpargs]
    elif isinstance(inpargs, dict):
        cmd = [gmx, backup_flag, "-nocopyright", gmxlib, *[i for tup in inpargs.items() for i in tup]]
    else:
        raise ValueError("inpargs must be either a list of input arguments or a dict")
    if interactive is None:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
    else:
        try:
            interactive = interactive.encode()
        except AttributeError:
            pass
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(interactive)
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    if proc.returncode != 0:
        print("Exited with error code", proc.returncode)
        if proc.returncode == 132:
            raise ChildProcessError("Core was dumped. This is probably due to an incompatible gromacs version")
        try:
            print(out.decode())
            print(err.decode())
        except UnboundLocalError:
            print("No error output")
        raise ChildProcessError(
            'Failed to execute command "{}" and input "{}"'\
            .format(' '.join(cmd), interactive))
    return out.decode(), err.decode()

def write_log(name_specifier, stdout, stderr, path="."):
    ''' Writes a logfile for gromacs output '''
    dt_str = datetime.now().strftime("%y.%m.%d-H%HM%MS%S")
    job_id = os.environ.get("SLURM_JOB_ID", None)

    if job_id is None:
        fname = "{}_{}.log".format(name_specifier, dt_str)
    else:
        fname = "slurmid{}_{}_{}.log".format(job_id, name_specifier, dt_str)
    with open(path + "/" + fname, "w") as flog:
        if job_id is not None:
            jname = os.environ["SLURM_JOB_NAME"]
            jnode = os.environ["SLURM_JOB_NODELIST"]
            flog.write("Job name: {}\nNodes: {}\n".format(jname, jnode))
        flog.write("STDOUT:\n")
        flog.write(stdout)
        flog.write("STDERR:\n")
        flog.write(stderr)



GMXNAME = find_executable("gmx_mpi")

