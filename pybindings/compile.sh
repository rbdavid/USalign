
#c++ -O3 -Wall -ffast-math -fvisibility=hidden -fPIC -shared -I "$(dirname "$PWD")" $(python3 -m pybind11 --includes) pySOIalign.cpp -o pySOIalign$(python3-config --extension-suffix)
c++ -O3 -Wall -fvisibility=hidden -fPIC -shared -I "$(dirname "$PWD")" $(python3 -m pybind11 --includes) pySOIalign.cpp -o pySOIalign$(python3-config --extension-suffix)
#c++ -O3 -Wall -shared -std=c++20 -fvisibility=hidden -I "$(dirname "$PWD")" -fPIC $(python3 -m pybind11 --includes) pySOIalign.cpp -o pySOIalign$(python3-config --extension-suffix)

#c++ -O3 -Wall -shared -std=c++20 -fvisibility=hidden -I "$(dirname "$PWD")" -fPIC $(python3 -m pybind11 --includes) pySOIalign.cpp -o pySOIalign$(python3-config --extension-suffix) -g

