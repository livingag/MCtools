import glob
import numpy as np
import pydicom
from scipy import interpolate
from matplotlib.path import Path
import os


def ctcreate(EGS_HOME, ctdir=".", structs={}):
    """
    Generates DOXSYZnrc .egsphant file using DICOM CT images and
    plan files.

    Interpolates CT images to the same dimensions and
    resolution as the planned dose grid. Sets any information outside
    of the `external` (Oncentra) and `body` (Eclipse) contours to
    0.
    :param EGS_HOME: location of EGS_HOME directory
    :param ctdir: directory containing DICOM CT and plan files
    :param structs: dictionary containing the name and densities of
                    any structures to override. e.g. {'Bolus': 1}
    :returns: .egsphant file in ctdir
    """

    files = []
    for fname in glob.glob(ctdir + "/RI*.dcm"):
        os.remove(fname)
    dosefile = glob.glob(ctdir + "/RD*.dcm")[0]
    rp = pydicom.dcmread(glob.glob(ctdir + "/RP*.dcm")[0])
    ID = pydicom.dcmread(dosefile).PatientID

    dose = DicomGrid(dosefile)
    dose.nz = len(dose.dicom.GridFrameOffsetVector)
    dose.dz = dose.dicom.GridFrameOffsetVector[1] - dose.dicom.GridFrameOffsetVector[0]

    zdim = (
        np.ones((1, dose.nz + 1)) * dose.z0
        - 0.5 * dose.dz
        + np.arange(0, dose.nz + 1) * dose.dz
    ) / 10
    zdim.sort()

    zct = (
        np.asarray(dose.dicom.GridFrameOffsetVector)
        + dose.dicom.ImagePositionPatient[2]
    )

    try:
        serial = pydicom.dcmread(glob.glob(ctdir + "/CT*.dcm")[0]).DeviceSerialNumber
        cttable = EGS_HOME + "/templates/{}.txt".format(serial)
        cttable = np.loadtxt(cttable)
    except:
        cttable = np.loadtxt(EGS_HOME + "/templates/ctcurve.txt")

    cted = interpolate.interp1d(cttable[:, 0], cttable[:, 1])
    edct = interpolate.interp1d(cttable[:, 1], cttable[:, 0])

    for filename in glob.glob(ctdir + "/CT*.dcm"):
        ct = pydicom.dcmread(filename)
        try:
            if ct.SliceLocation >= min(zct):
                if ct.PatientPosition == "FFS":
                    files.append((filename, -ct.SliceLocation))
                else:
                    files.append((filename, ct.SliceLocation))
        except:
            pass

    files.sort(key=lambda tup: tup[1])
    if not files:
        for z in zct:
            files.append((glob.glob(ctdir + "/CT*.dcm")[0], z))
    ct = DicomGrid(files[0][0])
    slices = np.asarray([i[1] for i in files])
    coords = [(x, y) for y in ct.yy for x in ct.xx]
    ctdim = int(np.sqrt(len(coords)))

    planct = np.zeros((len(files), ct.nx, ct.ny))

    try:
        rs = pydicom.dcmread(glob.glob(ctdir + "/RS*.dcm")[0])
        contours = {}
        for i in np.arange(len(rs.ROIContourSequence)):
            names = rs.StructureSetROISequence
            ROIname = names[i].ROIName
            contours[ROIname] = i
    except:
        pass

    couchmin = -np.inf
    couchxmin = np.inf
    couchxmax = -np.inf

    if all(x not in rp.Manufacturer for x in ["Varian"]):
        extid = [v for (k, v) in contours.items() if "external" in k.lower()]
        ext = rs.ROIContourSequence[extid[0]]
        try:
            couchid = [v for (k, v) in contours.items() if "couch" in k.lower()][0]
            couch = rs.ROIContourSequence[couchid]
            couchcont = np.array(couch.ContourSequence[0].ContourData).reshape((-1, 3))[
                :, :2
            ]

            cpath = Path(couchcont)
            cinpath = cpath.contains_points(coords).reshape((ctdim, ctdim))

            for sli in couch.ContourSequence:
                pts = np.array(sli.ContourData).reshape((-1, 3))
                if pts[:, 1].max() > couchmin:
                    couchmin = pts[:, 1].max()
                if pts[:, 0].min() < couchxmin:
                    couchxmin = pts[:, 0].min()
                if pts[:, 0].max() > couchxmax:
                    couchxmax = pts[:, 0].max()
        except:
            cinpath = np.full((ctdim, ctdim), False, dtype=bool)

    elif "Varian" in rp.Manufacturer:
        try:
            extid = [v for (k, v) in contours.items() if "external" in k.lower()][0]
        except:
            extid = [v for (k, v) in contours.items() if "body" in k.lower()][0]

        ext = rs.ROIContourSequence[extid]

        try:
            couchsurf = rs.ROIContourSequence[
                [v for (k, v) in contours.items() if "couchsurface" in k.lower()][0]
            ]
            couchint = rs.ROIContourSequence[
                [v for (k, v) in contours.items() if "couchinterior" in k.lower()][0]
            ]

            couchintcont = np.array(couchint.ContourSequence[0].ContourData).reshape(
                (-1, 3)
            )[:, :2]

            couchsurfcont = np.array(couchsurf.ContourSequence[0].ContourData).reshape(
                (-1, 3)
            )[:, :2]

            cintpath = Path(couchintcont)
            csurfpath = Path(couchsurfcont)

            cintinpath = cintpath.contains_points(coords).reshape((ctdim, ctdim))
            csurfinpath = csurfpath.contains_points(coords).reshape((ctdim, ctdim))

            for sli in couchsurf.ContourSequence:
                pts = np.array(sli.ContourData).reshape((-1, 3))
                if pts[:, 1].max() > couchmin:
                    couchmin = pts[:, 1].max()
                if pts[:, 0].min() < couchxmin:
                    couchxmin = pts[:, 0].min()
                if pts[:, 0].max() > couchxmax:
                    couchxmax = pts[:, 0].max()
        except:
            cintinpath = np.full((ctdim, ctdim), False, dtype=bool)
            csurfinpath = np.full((ctdim, ctdim), False, dtype=bool)

    phantom = False
    structconts = {}
    for (k, _) in contours.items():
        if any(
            x in k.lower() for x in ["bolus", "wire", "wax", "artifact", "artefact"]
        ):
            structs[k] = 1
    for key in structs.keys():
        structid = [v for (k, v) in contours.items() if k == key][0]
        struct = rs.ROIContourSequence[structid]
        structslices = {}
        for i, cont in enumerate(struct.ContourSequence):
            structslices.setdefault(cont.ContourData[2], []).append(i)
        structconts[key] = structslices

    contslices = {}
    for i, cont in enumerate(ext.ContourSequence):
        contslices.setdefault(cont.ContourData[2], []).append(i)

    for x, file in enumerate(files):
        ct = DicomGrid(file[0])
        ctz = file[1]
        ct.image = ct.image * ct.dicom.RescaleSlope + ct.dicom.RescaleIntercept
        inpath = np.full((ct.nx, ct.ny), False, dtype=bool)
        if ctz in contslices.keys():
            for contind in contslices[ctz]:
                contour = np.array(ext.ContourSequence[contind].ContourData).reshape(
                    (-1, 3)
                )[:, :2]
                cpath = Path(contour)
                inpath += cpath.contains_points(coords).reshape((ct.nx, ct.ny))
            ct.image[(inpath == True) & (ct.image < -953)] = -953
            ct.image[inpath == False] = -1000

        for key, structslices in structconts.items():
            structid = [v for (k, v) in contours.items() if k == key][0]
            if ctz in structslices.keys():
                inpath = np.full((ct.nx, ct.ny), False, dtype=bool)
                struct = rs.ROIContourSequence[structid]
                for structid in structslices[ctz]:
                    contour = np.array(
                        struct.ContourSequence[structid].ContourData
                    ).reshape((-1, 3))[:, :2]
                    cpath = Path(contour)
                    inpath += cpath.contains_points(coords).reshape((ct.nx, ct.ny))
                if not (inpath == False).all():
                    ct.image[inpath == True] = edct(structs[key])

        try:
            if all(x not in rp.Manufacturer for x in ["Varian", "Tomo"]):
                ct.image[cinpath == True] = 0
                planct[x, :, :] = ct.image
            elif "Varian" in rp.Manufacturer:
                ct.image[(cintinpath == False) & (csurfinpath == True)] = -300
                planct[x, :, :] = ct.image
        except:
            pass

    interpct = interpolate.interp1d(
        slices, planct, axis=0, bounds_error=False, fill_value=-1000
    )(zct)

    if (interpct[0, :, :] == -1000).all():
        x = 0
        while (interpct[x, :, :] == -1000).all():
            x += 1
        interpct[0:x, :, :] = interpct[x, :, :]

    if (interpct[-1, :, :] == -1000).all():
        x = -1
        while (interpct[x, :, :] == -1000).all():
            x -= 1
        interpct[x + 1 :, :, :] = interpct[x, :, :]

    densities = []
    materials = []

    if couchmin > dose.yy.max():
        rowadd = np.ceil((couchmin - dose.yy.max()) / dose.dy)
        dose.ny += rowadd
        dose.yy = np.append(dose.yy, np.arange(1, rowadd + 1) * dose.dy + dose.yy[-1])

    if couchxmin < dose.xx.min():
        coladd = np.ceil((dose.xx.min() - couchxmin) / dose.dx)
        dose.nx += coladd
        newx = np.arange(1, coladd + 1) * dose.dy + (
            dose.xx[0] - coladd * dose.dx - dose.dx
        )
        dose.xx = np.append(newx, dose.xx)
        dose.x0 = dose.xx[0]

    if couchxmax > dose.xx.max():
        coladd = np.ceil((couchxmax - dose.xx.max()) / dose.dx)
        dose.nx += coladd
        dose.xx = np.append(dose.xx, np.arange(1, coladd + 1) * dose.dx + dose.xx[-1])

    xdim = (
        np.ones((1, int(dose.nx) + 1)) * dose.x0
        - 0.5 * dose.dx
        + np.arange(0, int(dose.nx) + 1) * dose.dx
    ) / 10

    ydim = (
        np.ones((1, int(dose.ny) + 1)) * dose.y0
        - 0.5 * dose.dy
        + np.arange(0, int(dose.ny) + 1) * dose.dy
    ) / 10

    with open(ctdir + "/" + "%s_CT.egsphant" % (ID), "wb") as outfile:
        outfile.write(b"6\n")
        outfile.write(b"AIR521ICRU\n")
        outfile.write(b"LUNG521ICRU\n")
        outfile.write(b"ADIPOSE521ICRP\n")
        outfile.write(b"MUSCLE521ICRP\n")
        outfile.write(b"CARTILAGE521ICRU\n")
        outfile.write(b"ICRPBONE521ICRU\n")
        outfile.write(b"1 1 1 1 1 1\n")
        outfile.write(b"  %i  %i  %i\n" % (dose.nx, dose.ny, interpct.shape[0]))
        np.savetxt(outfile, np.sort(xdim), fmt="%.2f", delimiter=" ")
        np.savetxt(outfile, np.sort(ydim), fmt="%.2f", delimiter=" ")
        np.savetxt(outfile, zdim, fmt="%.2f", delimiter=" ")
        for x in np.arange(interpct.shape[0]):
            ct = DicomGrid(files[0][0])
            f = interpolate.interp2d(ct.xx, ct.yy, interpct[x, :, :])
            newct = f(
                dose.xx, dose.yy
            )  # * ct.dicom.RescaleSlope + ct.dicom.RescaleIntercept
            newct[newct > cttable[:, 0].max()] = cttable[:, 0].max()
            newcted = cted(newct)

            densities.append(cted(newct))

            for x in np.nditer(newcted, op_flags=["readwrite"]):
                if x < 0.0157:
                    x[...] = 1
                elif phantom:
                    x[...] = 2
                elif x < 0.58905:
                    x[...] = 2
                elif x < 0.98515:
                    x[...] = 3
                elif x < 1.07435:
                    x[...] = 4
                elif x < 1.35:
                    x[...] = 5
                else:
                    x[...] = 6

            materials.append(newcted)

        for x in materials:
            np.savetxt(outfile, x, fmt="%i", delimiter="")
            outfile.write(b"\n")

        for x in densities:
            for row in x:
                outfile.write(b" ")
                np.savetxt(outfile, row[np.newaxis], fmt="%.3f", delimiter=" ")
            outfile.write(b"\n")


class DicomGrid(object):
    def __init__(self, dicom):
        from pydicom import dcmread

        dicom = dcmread(dicom)
        self.dicom = dicom
        self.image = dicom.pixel_array
        self.x0 = dicom.ImagePositionPatient[0]
        self.dx = dicom.PixelSpacing[0]
        self.dy = dicom.PixelSpacing[1]
        self.nx = int(np.ceil(dicom.PixelSpacing[0] * dicom.Columns / self.dx))
        self.y0 = dicom.ImagePositionPatient[1]
        self.z0 = dicom.ImagePositionPatient[2]
        self.ny = int(np.ceil(dicom.PixelSpacing[1] * dicom.Rows / self.dx))
        self.dx *= dicom.ImageOrientationPatient[0]
        self.dy *= dicom.ImageOrientationPatient[4]
        self.xx = np.arange(self.x0, self.x0 + self.nx * self.dx, self.dx)
        self.yy = np.arange(self.y0, self.y0 + self.ny * self.dy, self.dy)
