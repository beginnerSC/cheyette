import subprocess, shutil, os

subprocess.check_call(['cmake', '-S', 'cheyette/src/_cheyette', '-B', 'build'])
subprocess.check_call(['cmake', '--build', 'build'])
shutil.copyfile('./build/Debug/_cheyette.cp312-win_amd64.pyd', './cheyette/_cheyette.pyd')
