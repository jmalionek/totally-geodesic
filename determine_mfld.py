import snappy, regina
import pandas as pd

def htlinkexterior12():
    with open('large_knots.txt', 'r') as data:
        knot12 = data.readlines()
    for i, knot in enumerate(knot12):
        knot = knot.strip('\n')
        if i < 6:
            knot = knot[0] + '_' + knot[1:]
        if 6 <= i < 40:
            knot = knot[:2] + '_' + knot[2:]
        knot12[i] = knot
    return knot12

def remove_bad_mfld(mfld_list):
    good_mfld = []
    for M_name in mfld_list:
        M = snappy.Manifold(M_name)
        if M.volume() > 12:
            try:
                M.link()
            except:
                good_mfld.append(M_name)
            else:
                if not M.link().is_alternating() and not M.is_two_bridge() :
                    good_mfld.append(M_name)
    return good_mfld

def htlinkexterior():
    chosen_mfld = []

    montesino = pd.read_csv('montesinos.csv', names=['name', 'tangle'])
    montesino_isosig = []
    for name in list(montesino['name']):
        M = snappy.Manifold(name)
        if M.volume() > 0.9:
            montesino_isosig.append(M.isometry_signature())

    for knot in htlinkexterior12():
        M_isosig = snappy.Manifold(knot).isometry_signature()
        if M_isosig not in montesino_isosig:
            chosen_mfld.append(knot)

    linkexterior = list(snappy.HTLinkExteriors(alternating=False, crossings=13)[11.9:])\
                   + list(snappy.HTLinkExteriors(alternating=False, crossings=14)[11.9:])\
                   + list(snappy.HTLinkExteriors(alternating=False, crossings=13)[11.9:])

    for link in linkexterior:
        try:
            link.isometry_signature()
        except:
            chosen_mfld.append(link.name())
        else:
            if link.isometry_signature() not in montesino_isosig:
                chosen_mfld.append(link.name())

    return chosen_mfld


if __name__ == '__main__':
    with open('htlinkexterior.txt', 'w') as f:
        mfldlist = htlinkexterior()
        for knot in mfldlist:
            f.write(knot + '\n')
