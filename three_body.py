#!/usr/bin/env python3
"""Simple three-body problem simulator.

This script simulates the motion of three bodies in a plane using
Newtonian gravity and a fourth-order Runge-Kutta integrator.

Example:
    python three_body.py --masses 5.972e24 5.972e24 7.348e22 \
        --positions -1e7,0 1e7,0 0,1.5e7 \
        --velocities 0,-1000 0,1000 800,0
"""
from __future__ import annotations

import argparse
from math import sqrt
from typing import Iterable, List, Sequence, Tuple

G = 6.67430e-11  # gravitational constant

Vector = Tuple[float, float]


# Predefined configurations that keep the bodies close enough for an
# interesting simulation. Selecting a preset overrides any manually
# specified masses, positions, and velocities.
PRESETS = {
    "demo": {
        "masses": [5.972e24, 5.972e24, 5.972e24],
        "positions": [(-1.0e7, 0.0), (1.0e7, 0.0), (0.0, 1.0e7)],
        "velocities": [(0.0, -1000.0), (0.0, 1000.0), (1000.0, 0.0)],
    },
    "figure8": {
        # Scaled version of the classic figure-eight three-body orbit.
        "masses": [5.972e24, 5.972e24, 5.972e24],
        "positions": [
            (-0.97000436e7, 0.24308753e7),
            (0.97000436e7, -0.24308753e7),
            (0.0, 0.0),
        ],
        "velocities": [
            (0.466203685e3, 0.43236573e3),
            (0.466203685e3, 0.43236573e3),
            (-0.93240737e3, -0.86473146e3),
        ],
    },
}


def v_add(a: Vector, b: Vector) -> Vector:
    return (a[0] + b[0], a[1] + b[1])


def v_sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1])


def v_mul(v: Vector, s: float) -> Vector:
    return (v[0] * s, v[1] * s)


def v_norm(v: Vector) -> float:
    return sqrt(v[0] * v[0] + v[1] * v[1])


def compute_accelerations(positions: Sequence[Vector], masses: Sequence[float]) -> List[Vector]:
    """Compute gravitational acceleration on each body."""
    acc: List[Vector] = []
    for i in range(3):
        ax, ay = 0.0, 0.0
        for j in range(3):
            if i == j:
                continue
            diff = v_sub(positions[j], positions[i])
            r = v_norm(diff)
            if r == 0:
                continue
            factor = G * masses[j] / (r ** 3)
            ax += diff[0] * factor
            ay += diff[1] * factor
        acc.append((ax, ay))
    return acc


def rk4_step(positions: Sequence[Vector], velocities: Sequence[Vector], masses: Sequence[float], dt: float) -> Tuple[List[Vector], List[Vector]]:
    """Perform a single Runge-Kutta step."""
    def advance(pos, vel, acc, h):
        return [v_add(pos[i], v_mul(vel[i], h)) for i in range(3)], \
               [v_add(vel[i], v_mul(acc[i], h)) for i in range(3)]

    a1 = compute_accelerations(positions, masses)
    p2, v2 = advance(positions, velocities, a1, dt / 2)
    a2 = compute_accelerations(p2, masses)
    p3, v3 = advance(positions, velocities, a2, dt / 2)
    a3 = compute_accelerations(p3, masses)
    p4, v4 = advance(positions, velocities, a3, dt)
    a4 = compute_accelerations(p4, masses)

    new_positions: List[Vector] = []
    new_velocities: List[Vector] = []
    for i in range(3):
        dp = v_mul(
            v_add(
                v_add(v_add(velocities[i], v_mul(v2[i], 2.0)), v_mul(v3[i], 2.0)),
                v4[i],
            ),
            dt / 6.0,
        )
        dv = v_mul(
            v_add(
                v_add(v_add(a1[i], v_mul(a2[i], 2.0)), v_mul(a3[i], 2.0)),
                a4[i],
            ),
            dt / 6.0,
        )
        new_positions.append(v_add(positions[i], dp))
        new_velocities.append(v_add(velocities[i], dv))
    return new_positions, new_velocities


def simulate(masses: Sequence[float], positions: Sequence[Vector], velocities: Sequence[Vector], dt: float, steps: int) -> List[List[Vector]]:
    """Simulate and return trajectory of each body."""
    traj: List[List[Vector]] = [[pos] for pos in positions]
    pos, vel = list(positions), list(velocities)
    for _ in range(steps):
        pos, vel = rk4_step(pos, vel, masses, dt)
        for i in range(3):
            traj[i].append(pos[i])
    return traj


def parse_vector(value: str) -> Vector:
    x_str, y_str = value.split(",")
    return (float(x_str), float(y_str))


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Three-body problem simulator")
    parser.add_argument(
        "--preset",
        choices=PRESETS.keys(),
        help="Use a predefined configuration with interesting initial conditions.",
    )
    parser.add_argument(
        "--masses",
        nargs=3,
        type=float,
        metavar=("M1", "M2", "M3"),
        default=[5.972e24, 5.972e24, 5.972e24],
        help="Masses of the bodies in kilograms.",
    )
    parser.add_argument(
        "--positions",
        nargs=3,
        type=parse_vector,
        metavar=("X1,Y1", "X2,Y2", "X3,Y3"),
        default=[(-1.0e7, 0.0), (1.0e7, 0.0), (0.0, 1.0e7)],
        help="Initial positions (x,y) in meters.",
    )
    parser.add_argument(
        "--velocities",
        nargs=3,
        type=parse_vector,
        metavar=("VX1,VY1", "VX2,VY2", "VX3,VY3"),
        default=[(0.0, -1000.0), (0.0, 1000.0), (1000.0, 0.0)],
        help="Initial velocities (vx,vy) in meters per second.",
    )
    parser.add_argument("--dt", type=float, default=1.0, help="Time step in seconds.")
    parser.add_argument("--steps", type=int, default=1000, help="Number of simulation steps.")
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.preset:
        preset = PRESETS[args.preset]
        args.masses = preset["masses"]
        args.positions = preset["positions"]
        args.velocities = preset["velocities"]

    trajectories = simulate(args.masses, args.positions, args.velocities, args.dt, args.steps)
    for i, body in enumerate(trajectories, start=1):
        x, y = body[-1]
        print(f"Body {i} final position: ({x:.3f}, {y:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
