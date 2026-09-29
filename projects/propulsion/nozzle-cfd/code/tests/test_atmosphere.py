"""US Standard Atmosphere 1976: layer-base values tabulated in the standard."""
import _path  # noqa: F401
from nozzlecfd.atmosphere import atmosphere, R_EARTH


def geometric(h):
    return R_EARTH * h / (R_EARTH - h)


def test_layer_bases():
    table = [  # geopotential m, T K, p Pa
        (0.0, 288.15, 101325.0),
        (11000.0, 216.65, 22632.06),
        (20000.0, 216.65, 5474.889),
        (32000.0, 228.65, 868.0187),
        (47000.0, 270.65, 110.9063),
    ]
    for h, T, p in table:
        Tc, pc, _ = atmosphere(geometric(h))
        assert abs(Tc - T) < 1e-6
        assert abs(pc - p) / p < 2e-5, (h, pc, p)


if __name__ == "__main__":
    test_layer_bases(); print("ok")
