import time
import numpy as np
import pandas as pd
from netCDF4 import Dataset

class pyton_netcdf(object):
    def __init__(self):
        pass


    def creat_netcdf_by_netcdf4(self,file_name_nc):
        ntimes = 5
        nlevels = 10
        nlats = 73
        nlons = 144


        rootgrp = Dataset(file_name_nc, "w", format="NETCDF4")
        print(rootgrp.data_model)

        rootgrp.description = "bogus example script"
        rootgrp.history     = "Created"+time.ctime(time.time())
        rootgrp.source      = "netDF4 pytohn module tutorial"

        time_dim  = rootgrp.createDimension(dimname="time",  size=None)
        level_dim = rootgrp.createDimension(dimname="level", size=None)
        lat_dim   = rootgrp.createDimension(dimname="lat",   size=nlats)
        lon_dim   = rootgrp.createDimension(dimname="lon",   size=nlons)


        times = rootgrp.createVariable(varname="time", datatype="f8", dimensions=("time",))
        times.units = "hours since 0001-01-01 00:00:00.0"
        times.calendar = "gregorian"

        levels = rootgrp.createVariable(varname="level", datatype="i4", dimensions=("level",))
        levels.units = "hPa"

        lats = rootgrp.createVariable(varname="lat", datatype="f8", dimensions=("lat",))
        lats.units = "degrees north"

        lons = rootgrp.createVariable(varname="lon", datatype="f8", dimensions=("lon",))
        lons.units = "degrees east"

        temp = rootgrp.createVariable(varname="temp", datatype="f8", dimensions=("time", "level", "lat", "lon",))
        #压缩
        #temp = rootgrp.createVariable(varname="temp", datatype="f8", ("temp",),zlib=True, least_significant_digit=3)
        temp.units ="mm"



        from datetime import datetime, timedelta
        from netCDF4 import date2num
        times_data = [datetime(2001,3,1)+n*timedelta(hours=12) for n in range(ntimes)]
        times[:]   = date2num(times_data, units=times.units, calendar = times.calendar)

        levels[:]  = np.arange(start=0, stop = nlevels, step=1)

        lats[:]    = np.arange(start=-90, stop=91, step=2.5)
        lons[:]    = np.arange(start=-180,stop=180,step=2.5)

        from numpy .random import uniform
        temp[:] = uniform(size=(ntimes,nlevels,nlats,nlons))

        rootgrp.close()
        print(file_name_nc)

    def read_netcdf_by_netcdf4(self,file_name_nc):
        rootgrp = Dataset(file_name_nc)

        #print(rootgrp.dimensions)
        print("dimensions")
        for dimobj in rootgrp.dimensions.values():
            #print(dimobj)
            if dimobj.isunlimited():
                print("\t\t", getattr(dimobj, "name"), "=UNLIMITED; //(", getattr(dimobj, "size"), "currently)")
                #也可以
                #print("\t\t", dimobj.name, "=UNLIMITED; //(", dimobj.size, "currently)")
            else:
                print("\t\t", dimobj.name, "=", dimobj.size)

        #print(rootgrp.variables)
        print("variables:")
        for varobj in rootgrp.variables.values():
            print("\t\t", varobj.dtype, varobj.name, varobj.dimensions, ":")
            for var_cattrs_name in varobj.ncattrs():
                print("\t\t"*2, varobj.name, ":", var_cattrs_name, "=", getattr(varobj,var_cattrs_name))

        print("//global attributes:")
        for global_attr in rootgrp.ncattrs():
            print("\t\t", global_attr, "=", getattr(rootgrp,global_attr))
        print("^"*20)

        # 直接读的例子
        print(rootgrp.description)
        print(rootgrp.dimensions["level"].name)
        print(rootgrp.dimensions["level"])
        print(rootgrp.variables["time"].units)
        print(rootgrp.variables["level"][:])
        print(rootgrp.variables["level"].shape)
        temp_nc_maskarray = rootgrp.variables["temp"][:]
        print(temp_nc_maskarray)
        print(temp_nc_maskarray.data)
        print(temp_nc_maskarray.mask)
        print(temp_nc_maskarray.fill_value)


    def netcdf_to_panel_by_netcdf4(self,file_name_nc):
        rootgrp = Dataset(file_name_nc)


        temp_panel = pd.Panel(data=rootgrp.variables["temp"][:, 0, :, :].data,
                              items= rootgrp.variables["time"][:].data,
                              major_axis=rootgrp.variables["lat"][:].data,
                              minor_axis=rootgrp.variables["lon"][:].data)
        print(temp_panel)

        from netCDF4 import num2date

        times = rootgrp.variables["time"]

        times_data = num2date(times[:], units=times.units, calendar=times.calendar)
        print(times_data)

        temp_panel4d = pd.Panel4D(data=rootgrp.variables["temp"][:, :, :].data,
                                labels=pd.period_range(start=times_data[0].strftime("%Y%m%d%H%M"), end=None, periods=times.shape[0],freq="12h",name="time"),
                                items=rootgrp.variables["level"][:].data,
                                major_axis=rootgrp.variables["lat"][:].data,
                                minor_axis=rootgrp.variables["lon"][:].data)
        print(temp_panel4d)










if __name__ == "__main__":

    py_nc = pyton_netcdf()
    #py_nc.creat_netcdf_by_netcdf4("netcdf4.nc")
    py_nc.read_netcdf_by_netcdf4("netcdf4.nc")
    py_nc.netcdf_to_panel_by_netcdf4("netcdf4.nc")






