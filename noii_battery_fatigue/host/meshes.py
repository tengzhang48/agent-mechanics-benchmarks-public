"""Quad4 meshes for the development problem D1 on [0,1]^2.

structured(n): n x n quads (n a multiple of 8), nodes counter-clockwise per element.
distorted(n, seed): the same mesh with interior nodes moved randomly by up to 0.25 h.
Nodes on the outer boundary and on the line y = 0.5 are not moved, so the notch
{0 <= x <= 0.125, y = 0.5} and the constraint points (1, 0.5), (1, 1) are nodes of every mesh.
"""
import numpy as np


def structured(n):
    if n % 8:
        raise ValueError("n must be a multiple of 8")
    xs = np.linspace(0.0, 1.0, n + 1)
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    nodes = np.column_stack([X.ravel(), Y.ravel()])
    elems = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            elems.append([a, a + 1, a + n + 2, a + n + 1])
    return nodes, np.asarray(elems, dtype=np.int64)


def distorted(n, seed=0, amplitude=0.25):
    nodes, elems = structured(n)
    h = 1.0 / n
    rng = np.random.default_rng(seed)
    x, y = nodes[:, 0], nodes[:, 1]
    fixed = (np.isclose(x, 0) | np.isclose(x, 1) | np.isclose(y, 0) | np.isclose(y, 1)
             | np.isclose(y, 0.5))
    move = ~fixed
    nodes = nodes.copy()
    nodes[move] += rng.uniform(-amplitude * h, amplitude * h, size=(move.sum(), 2))
    return nodes, elems


def notch_nodes(nodes, x_end=0.125):
    return np.where(np.isclose(nodes[:, 1], 0.5) & (nodes[:, 0] <= x_end + 1e-12))[0]


def boundary_nodes(nodes):
    x, y = nodes[:, 0], nodes[:, 1]
    return np.where(np.isclose(x, 0) | np.isclose(x, 1) | np.isclose(y, 0) | np.isclose(y, 1))[0]


def node_at(nodes, x, y):
    i = np.where(np.isclose(nodes[:, 0], x) & np.isclose(nodes[:, 1], y))[0]
    if len(i) != 1:
        raise ValueError(f"no unique node at ({x}, {y})")
    return int(i[0])
