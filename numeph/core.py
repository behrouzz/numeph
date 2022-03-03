import pickle, re
import numpy as np
from datetime import datetime
from jplephem import spk
from .julian import datetime_to_jd, jd_to_sec
from .utils import objects, num2txt, txt2num, str2seg


class SPK:
    """
    Load desired segments from bsp file to memory in desired interval
    
    Parameters
    ----------
        fname (str)     : path and name of bsp file
        t1 (datetime)   : ephemeris start time
        t2 (datetime)   : ephemeris end time
        segs_tup (list) : segments as (center, target) tuples
    """
    def __init__(self, fname, t1=None, t2=None, segs_tup=None):
        kernel = spk.SPK.open(fname)
        self.data = []

        # select segments
        all_segs = kernel.segments
        all_segs_tup = [(i.center, i.target) for i in all_segs]
        if segs_tup is None:
            segs_tup = all_segs_tup
            segments = all_segs
        else:
            segments = [i for i in all_segs if (i.center, i.target) in segs_tup]

        # select time interval
        slicing = False
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

            if coefficients.shape[1]==0:
                continue

            self.data.append([(seg.center, seg.target), domains, coefficients])
        kernel.close()

    def to_txt(self, fname):
        all_str = ''
        for seg_data in self.data:
            cen_tar, dom, coef = seg_data
            center, target = cen_tar

            cfx = coef[0,:,:]
            cfy = coef[1,:,:]
            cfz = coef[2,:,:]

            cf_xyz = np.concatenate((cfx,cfy,cfz))
            cf_str = num2txt(cf_xyz)

            n_cols = cf_xyz.shape[1]
            ini_dom = dom[0,0]
            dt_rec = dom[0,1] - dom[0,0]
            dt_dom = dom[-1,-1] - dom[0,0]
            n_recs = cfx.shape[0]

            first_row = [center, target, n_cols, n_recs, dt_rec, ini_dom, dt_dom]
            first_row = ['segment'] + [str(int(i)) for i in first_row]
            first_row = ','.join(first_row)+'\n'
            seg_str = first_row + cf_str
            all_str = all_str + seg_str
        f = open(fname, 'w')
        f.write(all_str)
        f.close()

    def to_pickle(self, fname):
        f = open(fname, 'wb')
        pickle.dump(self.data, f)
        f.close()


def load_txt(fname):
    """
    Load an ephemeris text file
    
    Parameters
    ----------
        fname (str) : path and name of the text file
    
    Returns
    ----------
        list: data as: [(center, target), domains, coefficients]
    """
    f = open(fname, 'r')
    all_str = f.read()
    f.close()
    fmt = 'segment,\d+,\d+,\d+,\d+,\d+,\d+,\d+\n'
    first_rows = re.findall(fmt, all_str)
    coef_rows = re.split(fmt, all_str)[1:]
    data = []
    for i in range(len(first_rows)):
        seg_data = str2seg(first_rows[i], coef_rows[i])
        data.append(seg_data)
    return data
    


def get_pos(file, seg_tup, t):
    """
    Get position of an object from a segment
    
    Parameters
    ----------
        file (str)      : path of pickle file
        seg_tup (tuple) : (center, target) of segment to be used
        t (datetime)    : time for which the position is requested
    
    Returns
    ----------
        pos (np.array): position of the object
    """
    
    f = open(file, 'rb')
    data = pickle.load(f)
    f.close()

    data = [i for i in data if i[0]==seg_tup][0]
    _, domains, coefficients = data
    
    jd = datetime_to_jd(t)
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


def geocentric(target, t, file):
    """
    Get geocentric position of an object
    
    Parameters
    ----------
        target (str) : name of the target, i.e. planet, sun or moon
        t (datetime) : time for which the position is requested
        file (str)   : path of pickle file
    
    Returns
    ----------
        pos (np.array): geocentric position of the object
    """
    
    target = target.lower()
    if target not in objects.keys():
        raise Exception('target not recognized!')
    target = objects[target]
    
    earB_ear = get_pos(file=file, seg_tup=(3,399), t=t)
    if target==301:
        earB_moo = get_pos(file=file, seg_tup=(3,301), t=t)
        pos = earB_moo - earB_ear
    elif target in [1,2,4,5,6,7,8,9,10]:
        SSB_plaB = get_pos(file=file, seg_tup=(0, target), t=t)
        SSB_earB = get_pos(file=file, seg_tup=(0,3), time=t)
        pos = SSB_plaB - earB_ear - SSB_earB
    return pos
