import numpy as np


def VMATmlcwrite(beam, filename, rmlc, mlcd):
    target = open(filename, 'w')
    target.write('%s\n' % filename)
    target.write('%s\n' % (beam.NumberOfControlPoints))

    for cp in beam.ControlPointSequence:
        try:
            mlcs = np.asarray(
                cp.BeamLimitingDevicePositionSequence[2].LeafJawPositions
            )*0.1
        except:
            mlcs = np.asarray(
                cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions
            )*0.1

        nmlc = beam[0x300a, 0xb6][2][0x300a, 0xbc].value

        for x, mlc in enumerate(mlcs):
            if mlc < 0:
                mlc = abs(mlc)
                neg = -1
            else:
                neg = 1
            mlcs[x] = (mlc*mlcd*0.01)+rmlc-(
                rmlc/mlcd * ((mlc*mlcd*0.01)**2 + mlcd**2)**0.5
            )
            mlcs[x] *= neg

        mlca = mlcs[:nmlc]
        mlcb = mlcs[nmlc:]

        mlcs = np.column_stack((mlca, mlcb))
        mu = float(cp.CumulativeMetersetWeight)

        mlcout = ''
        num = 1

        for x, mlc in reversed(list(enumerate(mlcs))):
            if x != 0 and (mlcs[x-1] == mlc).all():
                num += 1
            else:
                if num > 1:
                    mlcout = mlcout + '%.3f, %.3f, %i,\n' % (
                        mlc[0], mlc[1], num
                    )
                    num = 1
                else:
                    mlcout = mlcout + '%.3f, %.3f\n' % (mlc[0], mlc[1])
        target.write('%.5f\n' % mu)
        target.write('%s' % mlcout)

    target.close()


def VMATdynamlcwrite(filename, bankA, bankB):

    dlgA = np.loadtxt(bankA, skiprows=6, delimiter=',')
    dlgB = np.loadtxt(bankB, skiprows=6, delimiter=',')

    target = open(filename, 'w')
    target.write('%s\n' % filename)
    n = 0
    mu = 0
    for row in dlgA:
        if row[0] != mu:
            n += 1
            mu = row[0]
    target.write('%s\n' % (n))
    mu = 0

    for cpa, cpb in zip(dlgA, dlgB):
        if cpa[0] / 25000 != mu:
            mlcb = cpa[14:].reshape((-1, 4))[:, 1] * 0.001 + 0
            mlca = cpb[14:].reshape((-1, 4))[:, 1] * -0.001 + 0
            mlcs = np.column_stack((mlca, mlcb))
            mu = cpa[0] / 25000

            mlcout = ''
            num = 1

            for x, mlc in reversed(list(enumerate(mlcs))):
                if x != 0 and (mlcs[x-1] == mlc).all():
                    num += 1
                else:
                    if num > 1:
                        mlcout = mlcout + '%.3f,%.3f,%i\n' % (
                            mlc[0], mlc[1], num
                        )
                        num = 1
                    else:
                        mlcout = mlcout + '%.3f,%.3f\n' % (mlc[0], mlc[1])

            target.write('%.5f\n' % mu)
            target.write('%s' % mlcout)

    target.close()
