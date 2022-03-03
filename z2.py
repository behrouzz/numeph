f = '/home/behrouz/_py/Astronomy/Solar System/SPK/data/de440s.bsp'



from datetime import datetime
from numeph import SPK, load_txt

dc = load_txt('aaa.txt')

t = datetime.utcnow()

seg = dc[(3, 301)]

pos = seg.get_pos(t)
print(pos)
