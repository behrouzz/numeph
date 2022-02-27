**Author:** [Behrouz Safari](https://behrouzz.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

# numeph
*Convert JPL SPK ephemeris to numpy array*


## Installation

Install the latest version of *numeph* from [PyPI](https://pypi.org/project/numeph/):

    pip install numeph

Requirements are *numpy* and *jplephem*


## Save some segments from 'de440s.bsp' from 2020 to 2030:

Let's get the positions of the sun between two times:

```python
from numeph import save_segments

t1 = datetime(2020, 1, 1)
t2 = datetime(2030, 1, 1)

save_segments(in_file='de440s.bsp',
              out_file='de440s_2020_2030.pickle',
              t1=t1,
              t2=t2,
              seg_list=[*range(12)])
```

## get position of an object from a segment:

```python
from numeph import get_pos

t = datetime.utcnow()

pos = get_pos(file='de440s_2020_2030.pickle',
              segment_ind=2,
              time=t)
```

See more at [astrodatascience.net](https://astrodatascience.net/)
