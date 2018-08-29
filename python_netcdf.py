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
        '''
        temp_panel4d = pd.Panel4D(data=rootgrp.variables["temp"][:, :, :].data,
                                labels=pd.period_range(start=times_data[0].strftime("%Y%m%d%H%M"), end=None, periods=times.shape[0],freq="12h",name="time"),
                                items=rootgrp.variables["level"][:].data,
                                major_axis=rootgrp.variables["lat"][:].data,
                                minor_axis=rootgrp.variables["lon"][:].data)
        print(temp_panel4d)
        '''

    def creat_netcdf_by_xarray(self, file_name_nc):
        #Creat an example dateset
        import numpy as np
        import pandas as pd
        import xarray as xr

        ntimes = 5
        nlevels = 10
        nlats = 73
        nlons = 144

        levels  = np.arange(start=0, stop = nlevels, step=1)

        lats    = np.arange(start=-90, stop=91, step=2.5)
        lons    = np.arange(start=-180,stop=180,step=2.5)

        from numpy .random import uniform
        temps = uniform(size=(ntimes,nlevels,nlats,nlons))

        times = pd.date_range("2001-03-01 00:00", periods = temps.shape[0],freq="12h")

        dims = ["time", "level", "lat", "lon"]

        temp_attr = dict(standard_name="air_potential_tempeature", uints="f")
        prcp_attr = dict(standard_name="convective_precipitation_flux", uints="mm")

        level_attr = dict(standard_name="level", uints="hPa")
        lat_attr   = dict(standard_name="lat", uints="degrees north")
        lon_attr   = dict(standard_name="lon", uints="degrees east")

        ds = xr.Dataset({
            "temperature":(dims,temps,temp_attr),
            "precipitation":(dims,temps,prcp_attr)},
            coords={
                "time":times,
                "level":(["level"],levels,level_attr),
                "lat":(["lat"],lats,lat_attr),
                "lon":(["lon"],lons,lon_attr)
                #"reference_time":pd.Timestamp("2001-03-01 00:00")
                 }
            )

        ds.attrs=dict(description="bogus example script",
                    history="Created"+time.ctime(time.time()),
                    source="xarray python module tuorial",
                    author="some one")

        ds.to_netcdf(path=file_name_nc, mode="w", format="NETCDF4")

    def read_netcdf_to_xarray(self, file_name_nc):
        import xarray
        ds = xarray.open_dataset(file_name_nc)

        print("-"*20)

        #查看attrs，orderedDict，可以按照字典进行遍历
        print(ds.attrs)
        print("-"*20)

        #查看dims，SortedKeysDict，可以按照字典进行遍历
        print(ds.dims)
        print("-"*20)

        #查看coordinates 可以按照字典遍历，所谓coordinates，这里专指time level lat lon 这样的维度的数据
        print(ds.coords)
        print("-"*20)

        #coords_value 是一个DataArray
        for coords_value in ds.coords.values():
            print(coords_value)
            print("-"*5)
            print(coords_value.dims)
            print("-"*5)
            print(coords_value.values)
            print("-"*5)
            print(coords_value.coords)
            print("-"*5)
            print(coords_value.name)
            print("-"*5)
            print(coords_value.attrs)
            print("-"*5)
            print(coords_value.data)
        print("-"*20)

        #举个例子
        print(ds.coords["lon"].attrs)
        print(ds.coords["level"].values)

        print("-"*20)

        #查看数据，这个是除了维度信息后真正的数据，比如本例子中的温度和降雨量，也是个字典
        for data_value in ds.data_vars.values():
            print(data_value)
            print("-"*5)
            print(data_value.dims)
            print("-"*5)
            print(data_value.values)
            print("-"*5)
            print(data_value.coords)
            print("-"*5)
            print(data_value.name)
            print("-"*5)
            print(data_value.attrs)
            print("-"*5)
            print(data_value.data)

            #再举个例子
            print(ds.data_vars["precipitation"].values)  #同data
            print(ds.data_vars["precipitation"].attrs)   #参数
            print(ds.data_vars["precipitation"].name)    #名字
            print(ds.data_vars["precipitation"].data)    #真正的数据

            print(ds.data_vars["temperature"].values)
            print(ds.data_vars["temperature"].attrs)
            print(ds.data_vars["temperature"].name)
            print(ds.data_vars["temperature"].data)


        #这条命令类似于于ncdump -h
        print(ds.info())



if __name__ == "__main__":

    py_nc = pyton_netcdf()
    py_nc.creat_netcdf_by_netcdf4("netcdf4.nc")
    py_nc.read_netcdf_by_netcdf4("netcdf4.nc")
    py_nc.netcdf_to_panel_by_netcdf4("netcdf4.nc")

    py_nc.creat_netcdf_by_xarray("xarray.nc")
    py_nc.read_netcdf_to_xarray("xarray.nc")
