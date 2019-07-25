import pydicom
from pydicom.tag import Tag
import numpy as np
import re
import EDW, VMATmlc
from glob import glob
from numpy.random import random


def egsinpGenerate(plandir, tempdir, EGS_HOME, beamind=[]):

    BEAMGenerate(plandir, tempdir, EGS_HOME, beamind)
    DOSXYZGenerate(plandir, tempdir, EGS_HOME, beamind)


def BEAMGenerate(plandir, tempdir, EGS_HOME, beamind=[]):

    """
    Function for generation of input files for BEAMnrc.

    :param plandir: directory of where DICOM plan files are stored
    :param EGS_HOME: location of EGS_HOME directory
    :param beamind: list of indices of beams to generate files for - if
                    not specified, generate for all beams in plan
    :returns: BEAMnrc input files in plandir
    """

    planfile = glob(plandir + "/RP*.dcm")[0]

    plan = pydicom.dcmread(planfile)

    ID = plan.PatientID

    if not beamind:
        beamind = np.arange(len(plan.BeamSequence))

    for x in beamind:

        if (
            "kv" not in plan.BeamSequence[x].BeamName.lower()
            and "ELECTRON" not in plan.BeamSequence[x].RadiationType
        ):

            outfile = open(EGS_HOME + "/templates/BEAM_template.egsinp").read()

            xjawloc, yjawloc, rmlc, mlcd = getMachineInfo(outfile)

            if "DYNAMIC" in plan.BeamSequence[x].BeamType:
                btype = "Arc"

                beam = VMAT(plan.BeamSequence[x], xjawloc, yjawloc)

                fname = "%s_%s%i_BEAM.egsinp" % (ID, btype, beam.number)

                brem = (
                    max([abs(beam.yjaws[i]) for i in range(0, 2)]) ** 2
                    + max([abs(beam.xjaws[i]) for i in range(0, 2)]) ** 2
                ) ** 0.5 * 0.1 + 3

                mlcfile = "%s_%s%i.mlc" % (ID, btype, beam.number)

                if glob(plandir + "/*.dlg"):
                    bankA = sorted(glob(plandir + "/A*.dlg"))[x]
                    bankB = sorted(glob(plandir + "/B*.dlg"))[x]
                    VMATmlc.VMATdynamlcwrite(plandir + "/" + mlcfile, bankA, bankB)
                else:
                    VMATmlc.VMATmlcwrite(
                        plan.BeamSequence[x], plandir + "/" + mlcfile, rmlc, mlcd
                    )

                outfile = re.sub("!!W", "0", outfile)
                outfile = re.sub("!!X", beam.xout, outfile)
                outfile = re.sub("!!Y", beam.yout, outfile)

                outfile = re.sub("!!MLC", EGS_HOME + "/dosxyznrc/" + mlcfile, outfile)

                outfile = re.sub("!!modemlc", "1", outfile)

                outfile = re.sub("!!name", "%s Beam %i" % (ID, beam.number), outfile)
                outfile = re.sub("!!brem", "%.3f" % (brem), outfile)
                outfile = re.sub(
                    "!!rand", "%i, %i" % (random() * 31328, random() * 30081), outfile
                )

            else:

                beam = Beam(plan.BeamSequence[x], xjawloc, yjawloc, rmlc, mlcd)

                fname = "%s_Beam%i_BEAM.egsinp" % (ID, beam.number)

                brem = (
                    max([abs(beam.yjaws[i]) for i in range(0, 2)]) ** 2
                    + max([abs(beam.xjaws[i]) for i in range(0, 2)]) ** 2
                ) ** 0.5 * 0.1 + 3

                if beam.wedges > 0:
                    wedge = plan.BeamSequence[x].WedgeSequence[0]
                    filename = "%s_Beam%i_EDW.jaws" % (ID, beam.number)
                    edw = EDW.stt(
                        beam.energy,
                        wedge.WedgeAngle,
                        beam.yjaws[0] * 0.1,
                        beam.yjaws[1] * 0.1,
                        wedge.WedgeID[-1:],
                    )
                    EDW.EDWwrite(
                        edw,
                        plandir + "/" + filename,
                        beam,
                        wedge.WedgeID[-1:],
                        xjawloc,
                        yjawloc,
                    )

                    outfile = re.sub("!!W", "1", outfile)
                    outfile = re.sub("!!Y", "", outfile)
                    outfile = re.sub(
                        "" + str(yjawloc[0]) + ", " + str(yjawloc[1]) + ", \n",
                        "",
                        outfile,
                    )
                    outfile = re.sub(
                        "" + str(xjawloc[0]) + ", " + str(xjawloc[1]) + ", ",
                        "",
                        outfile,
                    )
                    outfile = re.sub(
                        "!!X", EGS_HOME + "/dosxyznrc/" + filename, outfile
                    )

                else:
                    outfile = re.sub("!!W", "0", outfile)
                    outfile = re.sub("!!X", beam.xout, outfile)
                    outfile = re.sub("!!Y", beam.yout, outfile)

                outfile = re.sub("!!name", "%s Beam %i" % (ID, beam.number), outfile)
                outfile = re.sub("!!brem", "%.3f" % (brem), outfile)
                outfile = re.sub(
                    "!!rand", "%i, %i" % (random() * 31328, random() * 30081), outfile
                )

                mlcout = ""
                num = 1

                for x, mlc in reversed(list(enumerate(beam.mlcs))):
                    if x != 0 and (beam.mlcs[x - 1] == mlc).all():
                        num += 1
                    else:
                        if num > 1:
                            mlcout = mlcout + "%.5f, %.5f, %i,\n" % (
                                mlc[0],
                                mlc[1],
                                num,
                            )
                            num = 1
                        else:
                            mlcout = mlcout + "%.5f, %.5f\n" % (mlc[0], mlc[1])

                outfile = re.sub("!!modemlc", "0", outfile)
                outfile = re.sub("!!MLC", mlcout[:-1], outfile)

            f_out = open(plandir + "/" + fname, "w")
            f_out.write(outfile)
            f_out.close()


def DOSXYZGenerate(plandir, tempdir, EGS_HOME, beamind=[]):

    """
    Function for generation of input files for DOSXYZnrc.

    :param plandir: directory of where DICOM plan files are stored
    :param EGS_HOME: location of EGS_HOME directory
    :param beamind: list of indices of beams to generate files for - if
                    not specified, generate for all beams in plan
    :returns: BEAMnrc input files in plandir
    """

    planfile = glob(plandir + "/RP*.dcm")[0]

    plan = pydicom.dcmread(planfile)

    ID = plan.PatientID

    if not beamind:
        beamind = np.arange(len(plan.BeamSequence))

    for x in beamind:

        if (
            "kv" not in plan.BeamSequence[x].BeamName.lower()
            and "ELECTRON" not in plan.BeamSequence[x].RadiationType
        ):

            if "DYNAMIC" in plan.BeamSequence[x].BeamType:

                outfile = open(
                    EGS_HOME + "/templates/DOSXYZ_template_VMAT.egsinp", "r"
                ).read()

                btype = "Arc"

                beam = VMAT(plan.BeamSequence[x])
                fname = "%s_%s%i_DOSXYZ.egsinp" % (ID, btype, beam.number)

                iso = beam.control[0].IsocenterPosition
                coll = 270 - beam.control[0].BeamLimitingDeviceAngle
                cpout = ""
                mu = 0

                if glob(plandir + "/*.dlg"):
                    bankA = sorted(glob(plandir + "/*.dlg"))[x]
                    dlg = np.loadtxt(bankA, skiprows=6, delimiter=",")

                    for cp in dlg:
                        if cp[0] / 25000 != mu:
                            cpout += "{0:.2f}, {1:.2f}, {2:.2f}, ".format(
                                iso[0] * 0.1, iso[1] * 0.1, iso[2] * 0.1
                            )
                            cpout += "90, "
                            cpout += "{0:.2f}, ".format(
                                (cp[6] / 10 + 180 - 2 * cp[6] / 10) - 90
                            )
                            cpout += "{0:.2f}, ".format(coll)
                            cpout += "45, "
                            mu = cp[0] / 25000
                            cpout += "{0:.5f}\n".format(cp[0] / 25000)
                else:
                    for cp in beam.control:
                        iso = beam.control[0].IsocenterPosition

                        if Tag(0x300A, 0x120) in cp.keys():
                            collimator = np.radians(cp[0x300A, 0x120].value)
                        if Tag(0x300A, 0x122) in cp.keys():
                            couch = np.radians(cp[0x300A, 0x122].value)
                        if Tag(0x300A, 0x011E) in cp.keys():
                            gantry = np.radians(cp[0x300A, 0x011E].value)

                        theta = np.degrees(np.arccos(-np.sin(couch) * np.sin(gantry)))
                        phi = np.degrees(
                            np.arctan2(-np.cos(gantry), np.cos(couch) * np.sin(gantry))
                        )
                        coll = np.degrees(
                            3 * np.pi / 2
                            - collimator
                            - np.arctan2(-np.sin(couch) * np.cos(gantry), np.cos(couch))
                        )
                        cpout += "{0:.4f}, {1:.4f}, {2:.4f}, ".format(
                            iso[0] * 0.1, iso[1] * 0.1, iso[2] * 0.1
                        )
                        cpout += "{0:.4f}, ".format(theta)
                        cpout += "{0:.4f}, ".format(phi)
                        cpout += "{0:.4f}, ".format(coll)
                        if btype == "Tomo":
                            cpout += "50, "
                        else:
                            cpout += "45, "
                        cpout += "{0:.10f}\n".format(cp.CumulativeMetersetWeight)

                outfile = re.sub(
                    "!!name", "{0} {1} {2}".format(ID, btype, beam.number), outfile
                )
                outfile = re.sub(
                    "!!CT", EGS_HOME + "/dosxyznrc/{}_CT.egsphant".format(ID), outfile
                )
                outfile = re.sub(
                    "!!BEAM", "{}_{}{}_BEAM".format(ID, btype, beam.number), outfile
                )
                outfile = re.sub(
                    "!!rand",
                    "{:.0f}, {:.0f}".format(random() * 31328, random() * 30081),
                    outfile,
                )
                outfile = re.sub("!!control", "{}".format(cpout[:-1]), outfile)

                if glob(plandir + "/*.dlg"):
                    n = 0
                    mu = 0
                    for row in dlg:
                        if row[0] != mu:
                            n += 1
                            mu = row[0]
                    outfile = re.sub("!!ncp", "{}".format(n), outfile)
                else:
                    outfile = re.sub(
                        "!!ncp",
                        "{}".format(plan.BeamSequence[x].NumberOfControlPoints),
                        outfile,
                    )

            else:

                beam = Beam(plan.BeamSequence[x])

                outfile = open(
                    EGS_HOME + "/templates/DOSXYZ_template.egsinp", "r"
                ).read()

                fname = "%s_Beam%i_DOSXYZ.egsinp" % (ID, beam.number)

                addcoll = 0
                addgantry = 0

                if plan.PatientSetupSequence[0].PatientPosition == "HFP":
                    addgantry = np.pi
                elif plan.PatientSetupSequence[0].PatientPosition == "FFS":
                    addcoll = np.pi

                outfile = re.sub(
                    "!!angle",
                    "%.2f"
                    % (
                        np.degrees(
                            np.arctan2(
                                -np.cos(beam.gantry + addgantry),
                                np.cos(beam.couch) * np.sin(beam.gantry + addgantry),
                            )
                        )
                    ),
                    outfile,
                )
                outfile = re.sub(
                    "!!coll",
                    "%.2f"
                    % (
                        np.degrees(
                            3 * np.pi / 2
                            - (beam.collimator + addcoll)
                            - np.arctan2(
                                -np.sin(beam.couch) * np.cos(beam.gantry + addgantry),
                                np.cos(beam.couch),
                            )
                        )
                    ),
                    outfile,
                )
                outfile = re.sub(
                    "!!couch",
                    "%.2f"
                    % (
                        np.degrees(
                            np.arccos(
                                -np.sin(beam.couch) * np.sin(beam.gantry + addgantry)
                            )
                        )
                    ),
                    outfile,
                )
                outfile = re.sub(
                    "!!iso",
                    "%.2f, %.2f, %.2f"
                    % (beam.iso[0] * 0.1, beam.iso[1] * 0.1, beam.iso[2] * 0.1),
                    outfile,
                )
                outfile = re.sub("!!name", "%s Beam %i" % (ID, beam.number), outfile)
                outfile = re.sub(
                    "!!CT", EGS_HOME + "/dosxyznrc/%s_CT.egsphant" % (ID), outfile
                )
                outfile = re.sub(
                    "!!BEAM", "%s_Beam%i_BEAM" % (ID, beam.number), outfile
                )
                outfile = re.sub(
                    "!!rand", "%i, %i" % (random() * 31328, random() * 30081), outfile
                )

            f_out = open(plandir + "/" + fname, "w")
            f_out.write(outfile)
            f_out.close()


def getBackscatter(plandir, template, beamind=[]):

    """
    Calculates backscatter correction for the beam specified due to
    backscatter from the upper jaws into the monitor chamber.

    :param plandir: directory of where DICOM plan files are stored
    :param beamind: index of beam
    :returns: backscatter factor to apply to dose calculations
    """

    plan = pydicom.dcmread(glob(plandir + "/RP*.dcm")[0])

    backscatter = {}

    if not beamind:
        beamind = np.arange(len(plan.BeamSequence))

    for x in beamind:

        if (
            "kv" not in plan.BeamSequence[x].BeamName.lower()
            and "ELECTRON" not in plan.BeamSequence[x].RadiationType
        ):

            xjawloc, yjawloc, rmlc, mlcd = getMachineInfo(template)

            if "DYNAMIC" in plan.BeamSequence[x].BeamType:
                beam = VMAT(plan.BeamSequence[x], xjawloc, yjawloc)
            else:
                beam = Beam(plan.BeamSequence[x], xjawloc, yjawloc, rmlc, mlcd)

            xjaws = [abs(x * 0.1) for x in beam.xjaws]
            yjaws = [abs(x * 0.1) for x in beam.yjaws]
            scatter = calcBackscatter(xjaws[0], xjaws[1], yjaws[0], yjaws[1])
            if beam.wedges == 0:
                backscatter[beam.number] = scatter
            else:
                backscatter[beam.number] = 1

    return backscatter


def calcBackscatter(x1, x2, y1, y2):

    R0 = 2.48287

    Rx = (0.8 - 0.0187 * (x1 + x2)) * (
        0.0395 * (y1 + y2) - 0.0000355 * (y1 ** 3 + y2 ** 3)
    )
    Ry = 3.08 - 0.0845 * (y1 + y2) + 0.0000447 * (y1 ** 3 + y2 ** 3)
    R = Rx + Ry

    Scb = (1 + 0.01 * R0) / (1 + 0.01 * R)

    return Scb


def getMachineInfo(template):

    """
    Gets static machine info from BEAMnrc template files

    :params energy: energy of beam
    :param plandir: location of plandir directory
    :returns xjawloc: top and bottom location of xjaw in model
    :returns yjawloc: top and bottom location of yjaw in model
    :returns rmlc: radius of mlc leaf in model
    :returns mlcd: centre of mlc leaf in model
    """
    # if energy == 6:
    #     fh = open(EGS_HOME+'/templates/BEAM_template.egsinp', 'r')

    # template = fh.read()
    # fh.close()
    yjawloc = re.search("Y\n(.*), !!Y", template).group(1).split(", ")
    yjawloc = [float(x) for x in yjawloc]
    xjawloc = re.search("X\n(.*), !!X", template).group(1).split(", ")
    xjawloc = [float(x) for x in xjawloc]
    mlczmin = float(re.search("MODE\n(.*), ZMIN\n(.*), ZTHICK", template).group(1))
    mlcthick = float(re.search("MODE\n(.*), ZMIN\n(.*), ZTHICK", template).group(2))
    mlcd = mlczmin + mlcthick / 2
    rmlc = float(re.search("\n(.*), ZFOCUS", template).group(1))

    return xjawloc, yjawloc, rmlc, mlcd


class Beam(object):
    def __init__(self, beam, xjawloc=None, yjawloc=None, rmlc=None, mlcd=None):
        control = beam[0x300A, 0x0111][0]
        self.number = beam.BeamNumber
        self.iso = control.IsocenterPosition
        self.energy = control.NominalBeamEnergy
        self.collimator = np.radians(control[0x300A, 0x120].value)
        self.couch = np.radians(control[0x300A, 0x122].value)
        self.wedges = beam.NumberOfWedges
        if xjawloc:
            self.xjaws = control[0x300A, 0x011A][0][0x300A, 0x011C].value
            self.x1 = [
                self.xjaws[0] * xjawloc[0] * 0.001,
                self.xjaws[0] * xjawloc[1] * 0.001,
            ]
            self.x2 = [
                self.xjaws[1] * xjawloc[0] * 0.001,
                self.xjaws[1] * xjawloc[1] * 0.001,
            ]
            self.xout = "%.3f, %.3f, %.3f, %.3f, " % (
                self.x2[0],
                self.x2[1],
                self.x1[0],
                self.x1[1],
            )
        if yjawloc:
            self.yjaws = control[0x300A, 0x011A][1][0x300A, 0x011C].value
            self.y2 = [
                self.yjaws[1] * -yjawloc[0] * 0.001,
                self.yjaws[1] * -yjawloc[1] * 0.001,
            ]
            self.y1 = [
                self.yjaws[0] * -yjawloc[0] * 0.001,
                self.yjaws[0] * -yjawloc[1] * 0.001,
            ]
            self.yout = "%.3f, %.3f, %.3f, %.3f, " % (
                self.y1[0],
                self.y1[1],
                self.y2[0],
                self.y2[1],
            )
        if rmlc:
            try:
                self.mlc = (
                    np.asarray(control[0x300A, 0x011A][2][0x300A, 0x011C].value) * 0.1
                )
                self.nmlc = beam[0x300A, 0xB6][2][0x300A, 0xBC].value
                for i, mlc in enumerate(self.mlc):
                    if mlc < 0:
                        mlc = abs(mlc)
                        neg = -1
                    else:
                        neg = 1
                    self.mlc[i] = neg * (
                        (mlc * mlcd * 0.01)
                        + rmlc
                        - (rmlc / mlcd * ((mlc * mlcd * 0.01) ** 2 + mlcd ** 2) ** 0.5)
                    )
                self.mlca = self.mlc[: self.nmlc]
                self.mlcb = self.mlc[self.nmlc :]
                self.mlcs = np.column_stack((self.mlca, self.mlcb))
            except:
                self.mlcs = np.column_stack((np.ones(60) * -20, np.ones(60) * 20))
        self.gantry = np.radians(control[0x300A, 0x011E].value)


class VMAT(object):
    def __init__(self, beam, xjawloc=None, yjawloc=None):
        ssd = beam.SourceAxisDistance
        self.control = beam[0x300A, 0x0111]
        self.number = beam.BeamNumber
        self.wedges = 0
        if xjawloc:
            self.xjaws = self.control[0][0x300A, 0x011A][0][0x300A, 0x011C].value
            self.x1 = [
                self.xjaws[0] * xjawloc[0] * 1 / ssd,
                self.xjaws[0] * xjawloc[1] * 1 / ssd,
            ]
            self.x2 = [
                self.xjaws[1] * xjawloc[0] * 1 / ssd,
                self.xjaws[1] * xjawloc[1] * 1 / ssd,
            ]
            self.xout = "%.4f, %.4f, %.4f, %.4f, " % (
                self.x2[0],
                self.x2[1],
                self.x1[0],
                self.x1[1],
            )
        if yjawloc:
            self.yjaws = self.control[0][0x300A, 0x011A][1][0x300A, 0x011C].value
            self.y2 = [
                self.yjaws[1] * -yjawloc[0] * 1 / ssd,
                self.yjaws[1] * -yjawloc[1] * 1 / ssd,
            ]
            self.y1 = [
                self.yjaws[0] * -yjawloc[0] * 1 / ssd,
                self.yjaws[0] * -yjawloc[1] * 1 / ssd,
            ]
            self.yout = "%.4f, %.4f, %.4f, %.4f, " % (
                self.y1[0],
                self.y1[1],
                self.y2[0],
                self.y2[1],
            )
