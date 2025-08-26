# 3-Body-Problem

Basic 3-body problem simulator.

## Running

Use `three_body.py` to simulate three planets on an x–y plane. The script
uses Newtonian gravity and lets you configure each planet's mass, initial
position, and velocity:

```
python three_body.py \
    --masses 5.972e24 5.972e24 7.348e22 \
    --positions -1e7,0 1e7,0 0,1.5e7 \
    --velocities 0,-1000 0,1000 800,0 \
    --dt 10 --steps 1000
```

All distances are in meters and velocities in meters per second. Adjust the
numbers to set different sizes (masses) and separations between the planets.
The program prints the final position of each planet after the simulation.

### Presets

Picking masses, positions, and velocities that produce an interesting
simulation can be tricky. The `--preset` option loads a predefined set of
values where the planets start close enough to interact visibly. For
example, to run the classic figure-eight choreography:

```
python three_body.py --preset figure8 --dt 1 --steps 2000
```

Presets override any manually specified `--masses`, `--positions`, and
`--velocities` arguments.
