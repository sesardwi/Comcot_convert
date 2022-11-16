import numpy as np
import time
import pygmt
import obspy

##---EDIT THIS AREA---
##---input file---
inp_path  = 'E:\\OTHER\\New_folder\\bu_weniza\\hasil_neartif\\i57j6mw8.7\\i53j10mw8.8\\'
hmax_path = inp_path +'hmax_layer05.dat'
lyrx_path = inp_path +'layer05_x.dat'
lyry_path = inp_path +'layer05_y.dat'
coast_path= inp_path +'coast_kupro.nc'
##---ouput file---
out_path  = 'E:\\OTHER\\New_folder\\bu_weniza\\hasil_neartif\\i57j6mw8.7\\i53j10mw8.8\\'
out_all_asc    = out_path +'i53j10mw8.8_all.asc'
out_all_xyz    = out_path +'i53j10mw8.8_all.xyz'
out_all_nc     = out_path +'i53j10mw8.8_all.nc'
out_land_asc   = out_path +'i53j10mw8.8_land.asc'
out_land_xyz   = out_path +'i53j10mw8.8_land.xyz'
out_land_nc    = out_path +'i53j10mw8.8_land.nc'
##---input parameter---
thr     = 0.025                             #inundation min. height threshold
spacing = "25e"                             #raster data spacing (e = meter)
##---------------------------------


t1=time.time()
##---read x- and y-layer---
layer_x = np.genfromtxt(lyrx_path)
layer_y = np.genfromtxt(lyry_path)
nx = len(layer_x)
ny = len(layer_y)
NN = nx*ny
xmin = np.min(layer_x)
ymin = np.min(layer_y)
xmax = np.max(layer_x)
ymax = np.max(layer_y)
dx = (layer_x[-1]-layer_x[1])/nx
dy = (layer_y[-1]-layer_y[1])/ny
cellsize = 0.5*(dx+dy)
region  = [xmin, xmax, ymin, ymax]      #output region area

##---read hmax file---
file = open(hmax_path,'r')
baris = file.readlines()
for i in range(len(baris)):
	baris[i]=baris[i].split()
file.close()

a_flt = np.ones(nx*ny)
no=0
for i in range(len(baris)):
    for j in range(len(baris[i])):
        if no < len(a_flt):
            a_flt[no]=float(baris[i][j])
            no=no+1

b_reshape = a_flt.reshape(ny,nx)
c_flip = np.flipud(b_reshape)

xyz_arr = np.zeros((NN,3))
xyz_data = b_reshape
for k in range(ny):
    for l in range(nx):
        q=(k)*nx
        xyz_arr[q+l][0] = layer_x[l]
        xyz_arr[q+l][1] = layer_y[k]
        xyz_arr[q+l][2] = xyz_data[k][l]

##---filter the data, get just on land and above inundation threshold---
inund_land     = pygmt.select(data=xyz_arr, region=region, mask='s/k',
                              gridmask=coast_path, z_subregion=str(thr)+'/100',resolution='h')
inund_land_arr = inund_land.to_numpy()


##---export hmax (land and sea) as nc file---
grd_inund_all  = pygmt.xyz2grd(data=xyz_arr, region=region, spacing=spacing,outgrid=out_all_nc)
##---export hmax (land) as nc file---
grd_inund_land = pygmt.xyz2grd(data=inund_land, region=region, spacing=spacing,outgrid=out_land_nc)


##---export hmax (land and sea) as xyz file---
outallxyz = open(out_all_xyz,'w')
for m in range(len(xyz_arr)):
    for n in range(len(xyz_arr[m])):
        if n == 0:
            outallxyz.write('%17.6f ' % xyz_arr[m][n])
        elif n == 1:
            outallxyz.write('%17.6f ' % xyz_arr[m][n])
        else:
            outallxyz.write('%8.3f\n' % xyz_arr[m][n])
outallxyz.close()

##---export hmax (land) as xyz file---
outlandxyz = open(out_land_xyz,'w')
for m in range(len(inund_land_arr)):
    for n in range(len(inund_land_arr[m])):
        if n == 0:
            outlandxyz.write('%17.6f ' % inund_land_arr[m][n])
        elif n == 1:
            outlandxyz.write('%17.6f ' % inund_land_arr[m][n])
        else:
            outlandxyz.write('%8.3f\n' % inund_land_arr[m][n])
outlandxyz.close()


##---export hmax (land and sea) as asc grid file---
outallasc = open(out_all_asc,'w')
outallasc.write('ncols %d \n' % (nx))
outallasc.write('nrows %d \n' % (ny))
outallasc.write('xllcorner %.10f \n' % (xmin))
outallasc.write('yllcorner %.10f \n' % (ymin))
outallasc.write('Cellsize %.10f \n' % (cellsize))
outallasc.write('Nodata_value -999.0000000000\n')
for m in range(len(c_flip)):
    for n in range(len(c_flip[m])):
        if c_flip[m][n] == 0:
            c_flip[m][n] = -999
        if n == 0:
            outallasc.write('%.8f' % c_flip[m][n])
        elif n == (nx-1):
            outallasc.write(' %.8f\n' % c_flip[m][n])
        else:
            outallasc.write(' %.8f' % c_flip[m][n])
outallasc.close()

##---export hmax (land) as asc grid file---
asc_xyz = pygmt.grd2xyz(grid=out_land_nc,nodata=-999, output_type='numpy',region=region)
minn = asc_xyz.min(axis=0)
maxx = asc_xyz.max(axis=0)
xfiltmin=minn[0]
yfiltmin=minn[1]
xfiltmax=maxx[0]
yfiltmax=maxx[1]
qq=np.zeros((len(asc_xyz),1))
for p in range(len(asc_xyz)):
    if asc_xyz[p][0]==asc_xyz[0][0]:
        qq[p]=p
    else:
        qq[p]=np.nan
qq=qq[~np.isnan(qq).all(axis=1)]
nnx=int(qq[1])
nny=int(len(asc_xyz)/nnx)
spacing2 = float(spacing[0:-1])
cllsz = obspy.geodetics.kilometer2degrees(spacing2/1000)

asc_z=np.zeros((len(asc_xyz),1))
for q in range(len(asc_xyz)):
    asc_z[q]=asc_xyz[q][2]
asc_rshp = asc_z.reshape(nny,nnx)

outlandasc = open(out_land_asc,'w')
outlandasc.write('ncols %d \n' % (nnx))
outlandasc.write('nrows %d \n' % (nny))
outlandasc.write('xllcorner %.10f \n' % (xfiltmin))
outlandasc.write('yllcorner %.10f \n' % (yfiltmin))
outlandasc.write('Cellsize %.10f \n' % (cllsz))
outlandasc.write('Nodata_value -999.000000\n')
for m in range(len(asc_rshp)):
    for n in range(len(asc_rshp[m])):
        if n == 0:
            outlandasc.write('%.8f' % asc_rshp[m][n])
        elif n == (nnx-1):
            outlandasc.write(' %.8f\n' % asc_rshp[m][n])
        else:
            outlandasc.write(' %.8f' % asc_rshp[m][n])
outlandasc.close()



t2=time.time() - t1
print('Elapsed time : '+str(t2)+' seconds')