"""."""
import os
import subprocess as sp

try:
    os.re('tbefile')
except OSError:
    pass
startupinfo = sp.STARTUPINFO()
startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
# Random output variable to avoid writing stuff from xfoil on the
# console
# sout = 0
# Calling xfoil with Popen
xfoil = sp.Popen(['xfoil.exe'],
                 stdin=sp.PIPE,
                 stdout=sp.PIPE,
                 stderr=sp.STDOUT,
                 startupinfo=startupinfo,
                 universal_newlines=True)

cmd = """
NORM
LOAD Rotor_root_bp3333
Rotor_root_bp3333
PANE
GDES
CADD




PANE
OPER
ITER 50
v
102920
MACH 0
PACC
tbefile

ALFA -29

QUIT
"""
# xfoil.stdin.write(cmd)
(out, err) = xfoil.communicate(cmd)
# for line in iter(xfoil.stdout.readline, b''):
#     sys.stdout.write(line)  # .decode(sys.stdout.encoding))
print(str(out))
print(cmd)
