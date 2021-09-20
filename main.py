from src.base_molecule import Cluster
from data.molecules_coordinates import (
    li,
    water,
    hydrogen2,
    hydrogen_fluoride,
    nano_boy,
    metal_complex,
)


def system_hf():
    hf = Cluster(hydrogen_fluoride)

    for i in range(20):
        ang = 45 * (i + 1)
        print(hf)
        hf = hf.rotate(0, x=0, y=0, z=45).translate(0, x=0.3, y=0, z=0)


def system_h2():
    h2 = Cluster(hydrogen2)

    for i in range(20):
        ang = 45 * (i + 1)
        print(h2)
        h2 = h2.rotate(0, x=0, y=0, z=45).translate(0, x=0.4, y=0, z=0)


def system_nano_boy():
    nb = Cluster(li, nano_boy)

    for i in range(20):
        angle = 90 * (-1) ** i
        print(nb)
        nb = nb.translate(1, x=0, y=0, z=0.5).rotate(1, x=angle, y=0, z=0)


def system_w_li():
    w = Cluster(li, water)

    for i in range(20):
        ang = 45 * (i + 1)
        print(w)
        w = w.translate(1, x=0.5, y=0, z=0).rotate(1, x=0, y=0, z=45)


def system_metal_complex():
    complex = Cluster(metal_complex)

    for i in range(20):
        print(complex)
        complex = complex.rotate(0, x=0, y=0, z=45).translate(0, x=2, y=0, z=0)


def system_w_li_cluster():
    from random import random, seed, randint

    w = Cluster(
        li,
        water,
        water,
        water,
        water,
        water,
        water,
    )

    total_steps = 1000
    seed(12345)
    for _ in range(total_steps):
        # molecule [0, w.total_fragments]
        mol = randint(0, w.total_fragments - 1)

        # angle between [0, 360)
        ax = random() * 360
        ay = random() * 360
        az = random() * 360

        # max_t: [-0.2, 0.2]
        max_t = 2
        tx = max_t * (random() - 0.5)
        ty = max_t * (random() - 0.5)
        tz = max_t * (random() - 0.5)

        print(w)
        w = w.translate(mol, x=tx, y=ty, z=tz).rotate(mol, x=ax, y=ay, z=az)


def test():
    w = Cluster(
        li,
        water,
    )

    li2 = Cluster(li).translate(0, 50, 0, 0)

    b = {
        "coordinates": [
            ("Li", 0.000000, 0.000000, 0.000000),
        ],
        "charge": +1,
        "multiplicity": 1,
    }

    w_li2 = Cluster(w, li2, b)

    print(w)
    print(li2)
    print(w_li2)


def run():
    # system_h2()
    # system_hf()
    # system_nano_boy()
    # system_w_li()
    # system_w_li_cluster()
    # system_metal_complex()

    # test
    test()


if __name__ == "__main__":
    run()
