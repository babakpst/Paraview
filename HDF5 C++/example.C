#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <hdf5.h>
 
// The number of cells in the X, Y dimensions
#define NX 30
#define NY 20
 

// ################################### Geometry
void
write_hdf5_geometry()
{
    hid_t     geometry_file_id;
    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[] = {NX+1, NY+1, 0};
    herr_t    status;
    const char *coordNames[] = {"/X", "/Y"};
 
    geometry_file_id = H5Fcreate("geometry.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
 
    // Create the coordinate data.
    float *x = (float *) malloc((NX+1)*(NY+1) * sizeof(float));
    float *y = (float *) malloc((NX+1)*(NY+1) * sizeof(float));
    int ndx = 0;
    for (int j = 0; j < NY+1; j++)
    {
        float yt = float(j) / float(NY);
        float angle = yt * M_PI;
        for (int i = 0; i < NX+1; i++)
        {
            float xt = float(i) / float(NX);
            float R = (1.-xt)*2. + xt*5.;
 
            x[ndx] = R * cos(angle);
            y[ndx] = R * sin(angle);
            ndx++;
        }
    }
 
    /* Write separate coordinate arrays for the x and y coordinates. */
    for(int did = 0; did < 2; ++did)
    {
        dims[0] = (NY + 1);
        dims[1] = (NX + 1);
        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate1(geometry_file_id, coordNames[did], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, did == 0 ? x : y);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
    }
 
    status = H5Fclose(geometry_file_id);

    free(x);
    free(y);
}



// ################################ Geometry
void
write_hdf5_timestep(int k)
{
    hid_t     data_file_id;
    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[] = {NX, NY, 0};
    herr_t    status;
    const char *coordNames[] = {"/X", "/Y"};

    char filename[256];
    sprintf(filename, "data-%d.h5", k);

    data_file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Create the scalar data.
    float *pressure = (float *) malloc(NX*NY * sizeof(float));
 
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int ndx = j * NX + i;
            pressure[ndx] = (float) (j + k);
        }
    }
 
    float *velocityx = (float *) malloc((NX+1)*(NY+1) * sizeof(float));
 
    for (int j = 0; j < NY+1; j++)
    {
        for (int i = 0; i < NX+1; i++)
        {
            int ndx = j * (NX+1) + i;
            velocityx[ndx] = (float) (i + k);
        }
    }
 
    // Write the scalar data.
    dims[0] = NY;
    dims[1] = NX;
    dataspace_id = H5Screate_simple(2, dims, NULL);
 
    dataset_id = H5Dcreate1(data_file_id, "/Pressure", H5T_NATIVE_FLOAT,
                           dataspace_id, H5P_DEFAULT);
 
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, pressure);
 
    status = H5Dclose(dataset_id);
 
    status = H5Sclose(dataspace_id);
 
    dims[0] = NY + 1;
    dims[1] = NX + 1;
    dataspace_id = H5Screate_simple(2, dims, NULL);
 
    dataset_id = H5Dcreate1(data_file_id, "/VelocityX", H5T_NATIVE_FLOAT,
                           dataspace_id, H5P_DEFAULT);
 
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, velocityx);
 
    status = H5Dclose(dataset_id);
 
    status = H5Sclose(dataspace_id);
 
    // Free the data.
    free(pressure);
    free(velocityx);
 
    status = H5Fclose(data_file_id);
}

void
write_xml_timestep(int k)
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh..
     */
    char filename[256];
    sprintf(filename, "xdmf_%06d.xmf", k);

    xmf = fopen(filename, "w");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "   <Time Value=\"%d\" />\n", k);
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        geometry.h5:/X\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        geometry.h5:/Y\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
    fprintf(xmf, "        data-%d.h5:/Pressure\n", k);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1);
    fprintf(xmf, "        data-%d.h5:/VelocityX\n", k);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fclose(xmf);
}

void
write_timestep(int k)
{
  write_hdf5_timestep(k);
  write_xml_timestep(k);
}
 

// ################################# Main 
int
main(int argc, char *argv[])
{
    write_hdf5_geometry();

    for (int i = 0; i < 100; i++)
      write_timestep(i);

    FILE *all = fopen("all.xmf", "w");
    fprintf(all, "<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n");
    fprintf(all, "  <Domain>\n");
    fprintf(all, "    <Grid GridType='Collection' CollectionType='Temporal'>\n");
    for (int i = 0; i < 100; i++)
      fprintf(all, "      <xi:include href='xdmf_%06d.xmf'/>\n", i);
    fprintf(all, "    </Grid>\n");
    fprintf(all, "  </Domain>\n");
    fprintf(all, "</Xdmf>\n");
      
    fclose(all);
   
    return 0;
}
