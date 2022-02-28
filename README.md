**Author:** [Behrouz Safari](https://behrouzz.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

# numeph
*Convert JPL SPK ephemeris to numpy array*


## Installation

Install the latest version of *numeph* from [PyPI](https://pypi.org/project/numeph/):

    pip install numeph

Requirements are *numpy* and *jplephem*


## Save some segments from 'de440s.bsp' from 2020 to 2030:

```python
from datetime import datetime
from numeph import save_segments

t1 = datetime(2020, 1, 1)
t2 = datetime(2030, 1, 1)

save_segments(in_file='de440s.bsp',
              out_file='de440s_2020_2030.pickle',
              t1=t1,
              t2=t2,
              segs_tup=[(0,10), (0,3), (3,399), (3,301)])
```

## get position of an object from a segment:

```python
from datetime import datetime
from numeph import get_pos

t = datetime.utcnow()

pos = get_pos(file='de440s_2020_2030.pickle',
              seg_tup=(0,3), # Earth Barycenter wrt SSB
              t=t)
```

## get geocentric position of an object:

```python
from datetime import datetime
from numeph import geocentric

t = datetime.utcnow()

pos = geocentric(target='moon', file='de440s_2020_2030.pickle', t=t)
```

See more at [astrodatascience.net](https://astrodatascience.net/)
