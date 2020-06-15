# SUMMARY:      create_tecplot.py
# USAGE:        create tecplot file for MODFLOW central plateau model
#               and extract parameters
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    December-2019
# DESCRIPTION:
# DESCRIP-END.
# COMMENTS:
#
# Last Change: 2020-03-10


import numpy as np
import re
from scipy.interpolate import RectBivariateSpline
import glob
import copy


def read_xy_plane(input_file):
    with open(input_file) as f:
        data = f.read()
    data = np.asarray([float(x) for x in re.split(" |\n", data) if x])
    # store z information
    data = data.reshape((nx, ny), order="F")
    data = data[:, ::-1]
    return(data)


def read_multi_refs(refs):
    multi_data = []
    for i_ref in refs:
        with open(i_ref) as f:
            data = f.read()
        data = np.asarray([float(x) for x in re.split(" |\n", data) if x])
        # store z information
        data = data.reshape((nx, ny), order="F")
        data = data[:, ::-1]
        multi_data.append(data)
    multi_data = np.swapaxes(np.swapaxes(
        np.asarray(multi_data), 0, 1), 1, 2)
    return(multi_data)


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


case_name = "P2Rv8.3"
case_name = "CPM_v845"

if case_name == "P2Rv8.3":
    case_dir = "/pic/projects/dvz/xhs_simus/Hanford_ACMs/" + \
        "Intera_PlateauToRiver_CP-57037-Rev2/calib_final/"
    nz = 7
    grid_file = case_dir+case_name+".spc"
    srf_refs = np.sort(glob.glob(case_dir+"/bot*")
                       )[::-1].tolist()+[case_dir+"top1.ref"]
    hy_refs = [case_dir+"kh" +
               str(iz+1)+".ref" for iz in range(nz)][::-1]  # 7 to 1
    kz_refs = [case_dir+"kh" +
               str(iz+1)+".ref" for iz in range(nz)][::-1]  # 7 to 1
    sy_refs = [case_dir+"sy" +
               str(iz+1)+".ref" for iz in range(nz)][::-1]  # 7 to 1
    sy_map = dict()
    kz_scale = 0.1
    output_dir = case_dir+"postprocess/"

elif case_name == "CPM_v845":
    case_dir = "/pic/projects/dvz/xhs_simus/Hanford_ACMs/" + \
        "Intera_CentralPlateau_CP-47631-Rev4/calibrated_model/historic_flow/"
    nz = 7
    grid_file = case_dir+"CPM.spc"
    srf_refs = [case_dir+"bot"+str(iz+1)+".ref" for iz in range(nz)][
        :: -1]+[case_dir+"top1.ref"]  # note the order is from 7 to 1
    hy_refs = [case_dir+"hy" +
               str(iz+1)+".ref" for iz in range(nz)][::-1]  # 7 to 1
    kz_refs = [case_dir+"hy" +
               str(iz+1)+".ref" for iz in range(nz)][::-1]  # 7 to 1
    sy_refs = [case_dir+"syzone_rng" +
               str(iz+1)+".inf" for iz in range(nz)][::-1]  # 7 to 1
    kz_scale = 1
    sy_map = dict()
    sy_map[1] = 2.5000000E-01
    sy_map[2] = 8.0036379E-02
    output_dir = case_dir+"postprocess/"
else:
    print("wrong case!")

# read grid information
with open(grid_file) as f:
    data = f.read()
data = [x for x in re.split("\n", data) if x]
ny, nx = map(int, data[0].split())
ox, ey = map(float, data[1].split()[0: 2])
data = [y for x in data[2:] for y in x.split()]
# convert string to grid spacing
delta = []
for s in data:
    if "*" in s:
        temp_n, temp_delta = s.split("*")
        delta += [float(temp_delta)]*int(temp_n)
    else:
        delta.append(float(s))
dx = np.array(delta[: nx])
dy = np.array(delta[-ny:])
oy = ey-np.sum(dy)
ex = ox+np.sum(dx)
x = ox+np.cumsum(dx)-0.5*dx
y = oy+np.cumsum(dy)-0.5*dy
cx = np.append(ox, ox+np.cumsum(dx))
cy = np.append(oy, oy+np.cumsum(dy))

# read surface
srf_data = read_multi_refs(srf_refs)
hy_data = read_multi_refs(hy_refs)
kz_data = read_multi_refs(kz_refs)
kz_data = kz_scale*kz_data
sy_data = read_multi_refs(sy_refs)
sy_ori = copy.deepcopy(sy_data)
for ikey in sy_map.keys():
    sy_data[sy_data == ikey] = sy_map[ikey]

# mark missing layer
z_array = np.zeros((nx, ny, nz))
for iz in range(nz):
    z_array[:, :, iz] = (srf_data[:, :, iz]+srf_data[:, :, iz+1])/2
    z_array[:, :, iz][srf_data[:, :, iz] == srf_data[:, :, iz+1]] = np.nan
    hy_data[:, :, iz][srf_data[:, :, iz] == srf_data[:, :, iz+1]] = -1
    kz_data[:, :, iz][srf_data[:, :, iz] == srf_data[:, :, iz+1]] = -1
    sy_data[:, :, iz][srf_data[:, :, iz] == srf_data[:, :, iz+1]] = -1


# interpolate surface to nodes (for use of Tecplot)
node_srf = []
for i_srf in range(srf_data.shape[-1]):
    f = RectBivariateSpline(x, y, srf_data[:, :, i_srf], kx=1, ky=1)
    node_srf.append(f(cx, cy))
node_srf = np.swapaxes(np.swapaxes(np.array(node_srf), 0, 1), 1, 2)

# output k and porosity
output_data = []
for iz in range(z_array.shape[2]):
    for iy in range(z_array.shape[1]):
        for ix in range(z_array.shape[0]):
            if ~np.isnan(z_array[ix, iy, iz]):
                output_data += [[x[ix], y[iy],
                                 z_array[ix, iy, iz],
                                 srf_data[ix, iy, iz],
                                 srf_data[ix, iy, iz+1],
                                 hy_data[ix, iy, iz],
                                 hy_data[ix, iy, iz],
                                 kz_data[ix, iy, iz],
                                 sy_data[ix, iy, iz]]]
output_data = np.array(output_data)

# output data
np.savetxt(output_dir+"model_grid_centroids.txt",
           output_data,
           fmt='%.8e',
           delimiter=',',
           header='Easting[m], Northing[m], Center Elevation [m, NAVD88],Bottom Elevation [m, NAVD88],Top Elevation [m, NAVD88], Ksxx[m/d], Ksyy[m/d], Kszz[m/d], Porosity')

# below is for tecplot
# prepare data for tecplot

varis = ["kz", "kh", "sy"]
varis_data = dict()
varis_data["kz"] = kz_data
varis_data["kh"] = hy_data
varis_data["sy"] = sy_data

# # variables
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
            cell_1 = west_ix+(south_iy-1)*(nx+1) + \
                (bottom_iz-1)*(nx+1)*(ny+1)
            # east_ix,south_iy,bottom_iz,
            cell_2 = east_ix+(south_iy-1)*(nx+1) + \
                (bottom_iz-1)*(nx+1)*(ny+1)
            # east_ix,north_iy,bottom_iz,
            cell_3 = east_ix+(north_iy-1)*(nx+1) + \
                (bottom_iz-1)*(nx+1)*(ny+1)
            # west_ix,north_iy,bottom_iz,
            cell_4 = west_ix+(north_iy-1)*(nx+1) + \
                (bottom_iz-1)*(nx+1)*(ny+1)
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

fname = open(output_dir+"tecplot.dat", 'w')
fname.write("\n".join(output_data))
fname.close()
