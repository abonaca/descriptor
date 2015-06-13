# descriptor
Package that creates observing scripts for the Dark Energy Camera

## Usage
Here we create a script in the local directory called 'M31_decam_script.json', containing three g band and five r band exposures of the Andromeda galaxy.

```python

import descriptor as desc
import astropy.coordinates as ac
import astropy.units as u

m31 = desc.Target('M31', ac.Angle('0:42:44 hours'), ac.Angle(41.269039*u.deg))

gband = desc.Sequence('M31_g', m31, exptype='object', exptime=270, band='g', 
                      seqtype='dither', ditherpattern='line', num=3)

rband = desc.Sequence('M31_r', m31, exptype='object', exptime=270, band='r', 
                      seqtype='dither', ditherpattern='center+rectangle', 
                      offset_ra=60*u.arcsec, offset_dec=90*u.arcsec)

	
script = desc.Script('M31_decam_script')
script.add_sequences([gband.exposures, rband.exposures])
script.write_json()
```
