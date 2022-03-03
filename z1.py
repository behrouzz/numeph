f = '/home/behrouz/_py/Astronomy/Solar System/SPK/data/de440s.bsp'



from datetime import datetime
from numeph import SPK

t1 = datetime(2020, 1, 1)
t2 = datetime(2030, 1, 1)

spk = SPK(fname=f, t1=t1, t2=t2,
          segs_tup=[(0,10), (0,3), (3,399), (3,301)])
"""
# save as txt file
spk.to_txt('data/de440s_2020_2030.txt')

# save as pickle
spk.to_pickle('data/de440s_2020_2030.pickle')
"""

t = datetime.utcnow()

seg = spk.segments[(3, 301)]

print(seg.get_pos(t))
spk.to_txt('aaa.txt')
