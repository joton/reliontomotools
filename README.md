# reliontomotools

Additional tools for subtomogram analysis in Relion tomo.


## Install
```bash
pip install reliontomotools
```

## Scripts

```python
warptomo2relion
```

Converts refinement of deformation models, particle poses and CTF paramaters from WARP/M to Relion tomo.

### Example
```bash
warptomo2relion -i 'WarpXML/TS_*.xml' -s 'tomograms/TS_*/TS_*_aligned.mrc' -d 1800 -o WarpConverted -p Refine3D/job010/run_data.star
```
