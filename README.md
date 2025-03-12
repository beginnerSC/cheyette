# PPC

Working example project putting together Poetry, Pybind11 and CMake. 

* One time set up: Run `poetry shell` and `poetry install`
* To build, run `poetry init` or `poetry build`
* Example usage:

```python
from ppc import add, sub

print(add(1, 2))
print(sub(1, 2))
```
* This project uses a `build.py` to trigger cmake run and when done copy the pyd file to the right place 
* A few things are hard coded which makes this project unportable:
    * In `./ppc/src/_ppc/CMakeLists.txt`: 
        * `set(pybind11_DIR C:/Users/chiuch3/AppData/Local/pypoetry/Cache/virtualenvs/ppc-MDAuejai-py3.12/Lib/site-packages/pybind11/share/cmake/pybind11)`
        * `include_directories(C:/Users/chiuch3/AppData/Local/pypoetry/Cache/virtualenvs/ppc-MDAuejai-py3.12/Lib/site-packages/pybind11/include)` to see pybind11
        * `include_directories(D:/UQL/python39-64/include)` to see Python.py
        * `link_directories(D:/UQL/python39-64/libs)` for linker to see python39.lib
    * In `./build.py`
        * `shutil.copyfile('./build/Debug/_ppc.cp312-win_amd64.pyd', './ppc/_ppc.pyd')`
* Where to change for new projects: 
    * all _ppc in `./ppc/src/_ppc/CMakeLists.txt`
    * `./ppc/src/_ppc/_ppc_pybind.cpp` which defines the module name in pybind
* F7 to build and F5 to debug c++ code. After code change always press F7 before F5 as `"preLaunchTask": "CMake: build"` in `launch.json` doesn't always work and is currently commented out 
    * If build files are gone, run `poetry build` or `CMake: Config` before pressing F7 to generate build files
* `pytest` to run pytest. Currently pytest builds the project in `setUpClass` so it takes 30 seconds