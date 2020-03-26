# SUMMARY:      create_tecplot.py
# USAGE:        create tecplot file for MODFLOW central plateau model
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    December-2019
# DESCRIPTION:
# DESCRIP-END.
# COMMENTS:
#
# Last Change: 2019-12-16


import numpy as np
import re
from scipy.interpolate import RectBivariateSpline


def wrap_line(values):
    """
    write long vectors in multiple line and limited digits
    """
    value_per_line = 10
    value_round = 6

    nline = int(np.ceil(len(values)/value_per_line))
    if nline == 1:
        return(["".join([str(np.round(s, value_round))+" " for s in values])])
    else:
        lines = []
        for iline in range(nline-1):
            line_value = values[: value_per_line]
            values = values[value_per_line:]
            lines += ["".join([str(np.round(s, value_round)) +
                               " " for s in line_value])]
        lines += ["".join([str(np.round(s, value_round)) +
                           " " for s in values])]
        return(lines)


# hard wired inputs start===================================
# file folder
simu_dir = "/pic/projects/dvz/HydroGeophys/plateau_model_fy20/"

# xy dimension
nx = 256
ny = 134
dx = np.full(nx, 100)
dy = np.full(ny, 100)
ox = 555600
oy = 129800
x = ox+np.cumsum(dx)-0.5*dx
y = oy+np.cumsum(dy)-0.5*dy
cx = np.append(ox, ox+np.cumsum(dx))
cy = np.append(oy, oy+np.cumsum(dy))

# z elevation files
nz = 7
srf_refs = ["top1.ref"]+["bot"+str(iz+1)+".ref" for iz in range(nz)]
srf_refs = srf_refs[::-1]

# variables
varis = ["kz", "hy", "kzone_rng", "syzone_rng"]
varis_extension = [".ref", ".ref", ".inf", ".inf"]

# hard wired inputs end===================================

# read surface
srf_data = []
for i_srf in srf_refs:
    ref_file = simu_dir+i_srf
    # read and reshape data
    with open(ref_file) as f:
        data = f.read()
    data = np.asarray([float(x) for x in re.split(" |\n", data) if x])
    # store z information
    data = data.reshape((nx, ny), order="F")
    data = data[:, ::-1]
    srf_data.append(data)

# interpolate surface to nodes
node_srf = []
for i_srf in srf_data:
    f = RectBivariateSpline(x, y, i_srf, kx=1, ky=1)
    node_srf.append(f(cx, cy))
node_srf = np.swapaxes(np.swapaxes(np.array(node_srf), 0, 1), 1, 2)

# read varis
varis_data = dict()
for ivari, iextension in zip(varis, varis_extension):
    varis_data[ivari] = []
    ivari_refs = [ivari+str(iz+1)+iextension for iz in range(nz)][::1]
    for i_ref in ivari_refs:
        ref_file = simu_dir+i_ref
        # read and reshape data
        with open(ref_file) as f:
            data = f.read()
        data = np.asarray([float(x) for x in re.split(" |\n", data) if x])
        # store vari data
        data = data.reshape((nx, ny), order="F")
        data = data[:, ::-1]
        varis_data[ivari].append(data)
    varis_data[ivari] = np.swapaxes(np.swapaxes(
        np.asarray(varis_data[ivari]), 0, 1), 1, 2)

output_data = ['TITLE = "FE-Volume Brick Data"']
output_data += ["VARIABLES = " +
                ",".join(['"'+x+'"' for x in ["X", "Y", "Z"] + varis])]
output_data += ['ZONE']
output_data += ['NODES='+str((nx+1)*(ny+1)*(nz+1))]
output_data += ['ELEMENTS='+str(nx*ny*nz)]
output_data += ['DATAPACKING=BLOCK']
output_data += ['ZONETYPE=FEBRICK']
output_data += ['VARLOCATION = ([1-3]= NODAL)']
output_data += ['VARLOCATION = ([4-'+str(len(varis)+7)+'] = CELLCENTERED)']
# put x, y, z in data
output_data += wrap_line(np.tile(cx, (ny+1)*(nz+1)))
output_data += wrap_line(np.tile(np.repeat(cy, nx+1), nz+1))
output_data += wrap_line(node_srf.flatten(order="F"))
# put varis in data
for ivari in varis:
    output_data += wrap_line(varis_data[ivari].flatten(order="F"))

# loop through cells
for iz in range(nz):
    bottom_iz = iz+1
    top_iz = iz+2
    for iy in range(ny):
        south_iy = iy+1
        north_iy = iy+2
        for ix in range(nx):
            west_ix = ix+1
            east_ix = ix+2
            # west_ix,south_iy,bottom_iz,
            cell_1 = west_ix+(south_iy-1)*(nx+1)+(bottom_iz-1)*(nx+1)*(ny+1)
            # east_ix,south_iy,bottom_iz,
            cell_2 = east_ix+(south_iy-1)*(nx+1)+(bottom_iz-1)*(nx+1)*(ny+1)
            # east_ix,north_iy,bottom_iz,
            cell_3 = east_ix+(north_iy-1)*(nx+1)+(bottom_iz-1)*(nx+1)*(ny+1)
            # west_ix,north_iy,bottom_iz,
            cell_4 = west_ix+(north_iy-1)*(nx+1)+(bottom_iz-1)*(nx+1)*(ny+1)
            # west_ix,south_iy,top_iz,
            cell_5 = west_ix+(south_iy-1)*(nx+1)+(top_iz-1)*(nx+1)*(ny+1)
            # east_ix,south_iy,top_iz,
            cell_6 = east_ix+(south_iy-1)*(nx+1)+(top_iz-1)*(nx+1)*(ny+1)
            # east_ix,north_iy,top_iz,
            cell_7 = east_ix+(north_iy-1)*(nx+1)+(top_iz-1)*(nx+1)*(ny+1)
            # west_ix,north_iy,top_iz,
            cell_8 = west_ix+(north_iy-1)*(nx+1)+(top_iz-1)*(nx+1)*(ny+1)
            output_data.append(" ".join(
                [str(x) for x in [cell_1, cell_2, cell_3, cell_4,
                                  cell_5, cell_6, cell_7, cell_8]]))

fname = open(simu_dir+"tecplot.dat", 'w')
fname.write("\n".join(output_data))
fname.close()
