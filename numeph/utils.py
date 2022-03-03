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


def str2seg(first_row, cf_str):
    first_row = first_row.split(',')[1:]
    first_row = [int(i) for i in first_row]
    center, target, n_cols, n_recs, dt_rec, ini_dom, dt_dom = first_row
    cf = txt2num(cf_str)
    cf = cf.reshape((3,n_recs,n_cols))
    domain = np.zeros((n_recs,2))
    domain[:,0] = [ini_dom + i*dt_rec for i in range(n_recs)]
    domain[:,1] = domain[:,0] + dt_rec
    seg_data = [(center,target), domain, cf]
    return seg_data
