import pickle
import numpy as np
from datetime import datetime
from jplephem.spk import SPK
from .julian import datetime_to_jd, jd_to_sec

objects = {'mercury':0, 'venus':1, 'mars':3, 'jupiter':4,
           'saturn':5, 'uranus':6, 'neptune':7, 'sun':8, 'moon':99}


def save_segments(in_file, out_file, t1=None, t2=None, seg_list=None):
    """
    Save segments of bsp file in desired interval as a pickle file
    
    Parameters
    ----------
        in_file (str)   : path and name of bsp file
        out_file (str)  : path and name of pickle file to be created
        t1 (datetime)   : ephemeris start time
        t2 (datetime)   : ephemeris end time
        seg_list (list) : index of segments to be included
    """
    
    kernel = SPK.open(in_file)
    data = []

    #select segments
    if seg_list is None:
        segments = kernel.segments
    else:
        segments = [kernel.segments[i] for i in seg_list]

    slicing = False
    
    # select time interval
    if (t1 is not None) and (t2 is not None):
        t1 = jd_to_sec(datetime_to_jd(t1))
        t2 = jd_to_sec(datetime_to_jd(t2))
        slicing = True
    
    for seg in segments:
        INIT, INTLEN, RSIZE, N = seg.daf.read_array(seg.end_i - 3, seg.end_i)
        t_ini, interval, coefficients = seg.load_array()

        cf_count = int(RSIZE - 2) // 3
        cf = seg.daf.map_array(seg.start_i, seg.end_i - 4)
        cf.shape = (int(N), int(RSIZE))

        MID_and_RADIUS = cf[:,:2]
        MID = MID_and_RADIUS[:,0]
        RADIUS = MID_and_RADIUS[:,1]

        domains = np.zeros(MID_and_RADIUS.shape)
        domains[:,0] = MID - RADIUS
        domains[:,1] = MID + RADIUS

        if slicing:
            mask1 = np.logical_and(t1>=domains[:,0], t1<domains[:,1])
            mask2 = np.logical_and(t2>=domains[:,0], t2<domains[:,1])
            rec1 = np.where(mask1)[0][0]
            rec2 = np.where(mask2)[0][0]
            coefficients = coefficients[:, rec1:rec2, :]
            domains = domains[rec1:rec2, :]

        data.append([domains, coefficients])

    kernel.close()
    f = open(out_file, 'wb')
    pickle.dump(data, f)
    f.close()



def get_pos(file, segment_ind, time):
    """
    Get position of an object from a segment
    
    Parameters
    ----------
        file (str): path of pickle file
        segment_ind (int): index of segment to be used
        time (datetime): time for which the position is requested
    
    Returns
    ----------
        pos (np.array): position of the object
    """
    
    f = open(file, 'rb')
    data = pickle.load(f)
    f.close()
    
    domains, coefficients = data[segment_ind]
    jd = datetime_to_jd(time)
    t = jd_to_sec(jd)
    
    mask = np.logical_and(t>=domains[:,0], t<domains[:,1])
    rec = np.where(mask)[0][0] # record index

    cfx = coefficients[0,rec,:]
    cfy = coefficients[1,rec,:]
    cfz = coefficients[2,rec,:]

    fx = np.polynomial.chebyshev.Chebyshev(coef=cfx, domain=domains[rec])
    fy = np.polynomial.chebyshev.Chebyshev(coef=cfy, domain=domains[rec])
    fz = np.polynomial.chebyshev.Chebyshev(coef=cfz, domain=domains[rec])

    pos = np.vstack((fx(t),fy(t),fz(t))).T[0]
    return pos


def geocentric(name, t, file):
    """
    Get geocentric position of an object
    
    Parameters
    ----------
        name (str)   : name of object (planets, sun or moon)
        t (datetime) : time for which the position is requested
        file (str)   : path of pickle file
    
    Returns
    ----------
        pos (np.array): geocentric position of the object
    """
    
    earB_ear = get_pos(file=file, segment_ind=10, time=t)
    if objects[name]==99:
        earB_moo = get_pos(file=file, segment_ind=9, time=t)
        pos = earB_moo - earB_ear
    else:
        SSB_plaB = get_pos(file=file, segment_ind=objects[name], time=t)
        SSB_earB = get_pos(file=file, segment_ind=2, time=t)
        pos = SSB_plaB - earB_ear - SSB_earB
    return pos
