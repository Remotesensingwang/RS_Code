import os.path
from osgeo import gdal_array
from osgeo import gdal
import numpy as np
from os import walk
import math
import os
import tarfile
'''
def unpacked(image_dir):
    f = []
    for (dirpath, dirnames, filenames) in walk(image_dir):
        f.extend(filenames)
        break
    for word in f:
       #if word.name.endswith('.tiff'):
        filename=image_dir+"\\"+word
        tf=tarfile.open(filename)
        tf.extractall('G:\data31')
        '''
class Landsat8Reader(object):
    def __init__(self,filename):
     self.bands = 7
     self.band_file_name = []
     self.nan_position = []
     self.filename=filename
     #self.file_name=os.path.join("G:\LAdata2019\lra",filename)

    def read(self):
     for band in range(self.bands):
         band_name =r'E:\IDLcode\python\LC08_L1GT_137032_20190819_20190902_01_T2_B1.TIF'
         self.band_file_name.append(band_name)

     ds = gdal.Open(self.band_file_name[0])
     image_dt = ds.GetRasterBand(1).DataType
     image = np.zeros((ds.RasterYSize, ds.RasterXSize, self.bands),
                      dtype=\
                       gdal_array.GDALTypeCodeToNumericTypeCode(image_dt))

     for band in range(self.bands):
         ds = gdal.Open(self.band_file_name[band])
         band_image = ds.GetRasterBand(1)
         a=band_image.ReadAsArray()
         image[:, :, band] = band_image.ReadAsArray()

     return image
     #print(image)

    def write(self, image, file_name, bands, format='Tiff'):
        ds = gdal.Open(self.band_file_name[0])
        projection = ds.GetProjection()
        geotransform = ds.GetGeoTransform()
        x_size = ds.RasterXSize
        y_size = ds.RasterYSize
        del ds
        band_count = bands
        dtype = gdal.GDT_Float32

        driver = gdal.GetDriverByName(format)
        new_ds = driver.Create(file_name, x_size, y_size, band_count, dtype)
        new_ds.SetGeoTransform(geotransform)
        new_ds.SetProjection(projection)

        for band in range(self.bands):
            new_ds.GetRasterBand(band + 1).WriteArray(image[:, :, band])
            new_ds.FlushCache()
        del new_ds

    def radiometric_calibration(self):
        image = self.read()

        def get_calibration_parameters():
            # filename1=self.filename+ "\\"+os.path.basename(filename)+ "_MTL" + ".txt"
            filename1=r'E:\IDLcode\python\LC08_L1GT_137032_20190819_20190902_01_T2_MTL.txt'
            f = open(filename1, 'r')
            metadata = f.readlines()
            f.close()
            multi_parameters = []
            add_parameters = []
            sun_e=0
            parameter_start_line = 0

            for lines in metadata:
                test_line = lines.split('=')
                if test_line[0] == 'RADIANCE_MULT_BAND_1 ':
                    break
                else:
                    parameter_start_line = parameter_start_line + 1

            for lines in range(parameter_start_line-1, parameter_start_line + 11):
                parameter = float(metadata[lines].split('=')[1])
                multi_parameters.append(parameter)

            for lines in range(parameter_start_line + 11, parameter_start_line + 22):
                parameter = float(metadata[lines].split('=')[1])
                add_parameters.append(parameter)

            for lines in metadata:
                test_line = lines.split('=')
                if test_line[0] == 'SUN_ELEVATION':
                    sun_e=float(test_line[1])
            #print(multi_parameters, add_parameters)
            return multi_parameters, add_parameters, sun_e

        multi_parameters, add_parameters,sun_e = get_calibration_parameters()
        cali_image = np.zeros_like(image,dtype=float)
        toa=np.zeros_like(image,dtype=float)

        for band in range(self.bands):
            gain = multi_parameters[band]
            offset = add_parameters[band]
            cali_image[:, :, band] = image[:, :, band] * gain + offset
            toa = cali_image[:, :, band]/math.sin(sun_e)
        del image
        return toa


if __name__ == "__main__":
    # gdal.AllRegister()
    #image_dir = r"G:\data3"
    #unpacked(image_dir)

    image_dir=r"E:\IDLcode\python\data"
    f = []
    for (dirpath, dirnames, filenames) in walk(image_dir):
        f.extend(filenames)
        break
    for word in f:
        filename = image_dir + "\\" + word
        data = Landsat8Reader(filename)
        image = data.read()
        toa = data.radiometric_calibration()
        file_path = r'E:\data'
        data.write(toa, file_path, data.bands)
