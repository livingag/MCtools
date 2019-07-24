import numpy as np


class stt(object):

    def __init__(self, energy, angle, y1, y2, jaw):
        self.energy = energy
        self.angle = angle
        self.y1 = y1
        self.y2 = y2
        self.position = np.linspace(-20, 10, 31)  # Golden STT table positions
        if energy == 6:
            self.dose = np.array([0.150691, 0.168051, 0.187220, 0.208376,
                                  0.231707, 0.257422, 0.285748, 0.316933,
                                  0.351245, 0.388978, 0.430453, 0.476017,
                                  0.526050, 0.580964, 0.641210, 0.707274,
                                  0.779690, 0.859035, 0.945937, 1.041080,
                                  1.145206, 1.259122, 1.383704, 1.519904,
                                  1.668756, 1.831381, 2.008999, 2.202931,
                                  2.414611, 2.645597, 2.897577])

            # Backscatter correction factors
            a = 1.03
            b = 7.3e-4
        else:
            self.dose = np.array([0.195441, 0.214977, 0.236228, 0.259328,
                                  0.284425, 0.311673, 0.341242, 0.373311,
                                  0.408074, 0.445738, 0.486525, 0.530673,
                                  0.578438, 0.630093, 0.685931, 0.746266,
                                  0.811433, 0.881793, 0.957731, 1.039658,
                                  1.128016, 1.223276, 1.325944, 1.436558,
                                  1.555697, 1.683978, 1.822059, 1.970646,
                                  2.130494, 2.302406, 2.487244])

            # Backscatter correction factors
            a = 1.02
            b = 4.9e-4

        # STT calculation from Golden
        angle = np.radians(angle)
        w0 = (np.tan(np.radians(60))-np.tan(angle))/(np.tan(np.radians(60)))
        w60 = np.tan(angle)/np.tan(np.radians(60))
        self.dose = self.dose*w60 + self.dose[20]*w0

        # Backscatter correction
        self.dose = (a-b*np.abs(self.position))*self.dose

        # EDW jaw positions (every EDW STT has 20 jaw positions)
        if jaw == '1':
            yend = y2-0.5
            ypos = np.linspace(y1, y2-0.5, 20)
            i = 1
        else:
            yend = y1+0.5
            ypos = -np.linspace(y2, y1+0.5, 20)
            i = -1

        # STT interpolation / normalisation
        self.dose = np.interp(ypos, self.position, self.dose)
        self.position = i*ypos
        self.dose = (self.dose / self.dose.max())


def EDWwrite(stt, filename, beam, jaw, xjawloc, yjawloc):
    # Input file writing
    target = open(filename, 'w')
    target.write('%s MV %s deg wedge\n' % (stt.energy, stt.angle))
    target.write('20\n')
    for i in range(0, 20):
        target.write('%.4f\n' % stt.dose[i])
        if jaw == '1':
            y1t = -(28.0/100)*stt.position[i]
            y1b = -(35.8/100)*stt.position[i]
            y2t = beam.y2[0]
            y2b = beam.y2[1]
        else:
            y2t = -(28.0/100)*stt.position[i]
            y2b = -(35.8/100)*stt.position[i]
            y1t = beam.y1[0]
            y1b = beam.y1[1]
        target.write('%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n' % (28.0, 35.8, y1t, y1b, y2t, y2b))
        target.write('%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n' % (36.7, 44.5, beam.x2[0], beam.x2[1], beam.x1[0], beam.x1[1]))
    target.close()
