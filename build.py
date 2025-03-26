import subprocess, shutil

subprocess.check_call(['cmake', '-S', 'cheyette/src/_cheyette', '-B', 'build'])
subprocess.check_call(['cmake', '--build', 'build'])
shutil.copyfile('./build/Debug/_cheyette.cp312-win_amd64.pyd', './cheyette/_cheyette.pyd')

# subprocess.check_call(['poetry', 'run', 'pip', 'uninstall', '-y', 'cheyette'])
# subprocess.check_call(['poetry', 'run', 'pip', 'install', './dist/cheyette-0.1.0-cp312-cp312-win_amd64.whl'])     # with these two lines, wheel will fail to build (of course)