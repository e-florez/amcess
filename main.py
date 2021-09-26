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

    print(hf)
    print(hf.total_atoms)
    print(hf.total_fragments)
    print(hf.symbols)
    print(hf.total_mass)
    print(hf.center_of_mass)
    print(hf.principal_axis)
    print(hf.add_fragments(water))
    print(hf.delete_fragments(4))

    w = hf.add_fragments(water)
    print(w.translate(1, x=20, y=10, z=50))
    print(w.rotate(1, x=0, y=0, z=45))
    print()

    # for i in range(20):
    #     ang = 45 * (i + 1)
    #     print(hf)
    #     hf = hf.rotate(0, x=0, y=0, z=45).translate(0, x=0.3, y=0, z=0)


def system_h2():
    h2 = Cluster(hydrogen2)

    for i in range(20):
        print(h2)
        h2 = h2.rotate(0, x=0, y=0, z=45).translate(0, x=0.4, y=0, z=0)


def system_nano_boy():
    nb = Cluster(li, nano_boy)

    with open("nano_boy.xyz", "w") as file_xyz:
        #     file_xyz.write(str(nb))

        for i in range(2):
            angle = 90 * (-1) ** i
            print(nb)
            nb = nb.translate(1, x=0, y=0, z=0.5).rotate(1, x=angle, y=0, z=0)

            # file_xyz.write(str(nb))


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
    # from random import random, seed, randint
    from numpy import random

    w = Cluster(
        li,
        water,
        water,
        water,
        water,
        water,
        water,
    )

    # print(w)

    total_steps = 1000

    random_gen = random.default_rng(1234)

    probability = [
        1 / (w.total_fragments - 1) if p else 0 for p in range(w.total_fragments)
    ]

    with open("w6_li.xyz", "w") as file_xyz:
        for _ in range(total_steps):
            # molecule [0, w.total_fragments]
            mol = random_gen.choice(
                w.total_fragments,
                p=probability,
                # p=[0, 1 / 6, 1 / 6, 1 / 6, 1 / 6, 1 / 6, 1 / 6],
            )

            # angle between [0, 360)
            ax = random_gen.uniform() * 360
            ay = random_gen.uniform() * 360
            az = random_gen.uniform() * 360

            # max_t: [-2, 2]
            max_t = 2
            tx = max_t * (random_gen.uniform() - 0.5)
            ty = max_t * (random_gen.uniform() - 0.5)
            tz = max_t * (random_gen.uniform() - 0.5)

            # print(w)
            w = w.translate(mol, x=tx, y=ty, z=tz).rotate(mol, x=ax, y=ay, z=az)

            file_xyz.write(w.xyz)


def test():
    w = Cluster(
        li,
        water,
    )

    li2 = Cluster(li).translate(0, 50, 0, 0)

    b = {
        "atoms": [
            ("Xe", 0.000000, 0.000000, 0.000000),
        ],
        "charge": +1,
        "multiplicity": 1,
    }

    # w_li2 = Cluster(w, li2, b)

    # with open("abc.xyz", "w") as file_xyz:
    #     file_xyz.write(str(w))

    print(w)
    print("str" + "-" * 10)
    print(str(w))
    # str(w)
    # print("repr" + "-" * 10)
    # print(repr(w))
    # # print(li2)
    # # print(w_li2)


def run():
    pass
    # system_h2()
    system_hf()
    # system_nano_boy()
    # system_w_li()
    # system_w_li_cluster()
    # system_metal_complex()

    # test
    test()


if __name__ == "__main__":
    run()
