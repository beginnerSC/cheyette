import subprocess, shutil, os

subprocess.check_call(['cmake', '-S', 'ppc/src/_ppc', '-B', 'build'])
subprocess.check_call(['cmake', '--build', 'build'])
shutil.copyfile('./build/Debug/_ppc.cp312-win_amd64.pyd', './ppc/_ppc.pyd')
