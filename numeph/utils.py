import io
import numpy as np

objects = {'solar system barycenter': 0,
           'mercury barycenter': 1,
           'venus barycenter': 2,
           'earth barycenter': 3,
           'mars barycenter': 4,
           'jupiter barycenter': 5,
           'saturn barycenter': 6,
           'uranus barycenter': 7,
           'neptune barycenter': 8,
           'pluto barycenter': 9,
           'sun': 10,
           'moon': 301,
           'earth': 399,
           'mercury': 199,
           'venus': 299}


def num2txt(arr):
    output = io.BytesIO()
    np.savetxt(output, arr)
    return output.getvalue().decode('utf-8')

def txt2num(arr_str):
    arr_b = arr_str.encode('utf-8')
    return np.loadtxt(io.BytesIO(arr_b))
