#include <string.h>
#include <iostream>
#include <fstream>

#include <Xdmf.h>
#include <XdmfDOM.h>
#include <XdmfDomain.h>
#include <math.h>

#undef STRUCTURED
#define HEAVY 1

#define NX 8
#define NY 8
#define NZ 8

int np = 2;

int
doit(const char *name, int dir)
{
  char dsn[128];
  int psz = NX / np;

  XdmfDOM       dom = XdmfDOM();
  XdmfRoot      root = XdmfRoot();

  root.SetDOM(&dom);
  root.SetVersion(2.2);
  root.Build();

  XdmfDomain    domain = XdmfDomain();
  root.Insert(&domain);

  XdmfGrid *pgrid = new XdmfGrid;
  pgrid->SetGridTypeFromString("Collection");
  pgrid->SetCollectionTypeFromString("Spatial");
  domain.Insert(pgrid);

  for (int p = 0; p < np; p = p + 1)
  {
    int sx = p * psz;
    int ex = (p == np-1) ? NX-1 : (p+1) * psz;
    int nx = (ex - sx) + 1;
    XdmfInt32 indx = 0;

    XdmfGrid  *grid = new XdmfGrid;
    pgrid->Insert(grid);

    XdmfTopology *top = grid->GetTopology();
    XdmfGeometry *geom = grid->GetGeometry();

#ifdef STRUCTURED
    XdmfInt32 nd = 3;
    XdmfInt64 cdim[] = {(nx-1), (NY-1), (NZ-1)};

    double orig[] = {-1.0 + (2.0*sx/(NX-1)), -1.0, -1.0};
    double deltas[] = {2.0/(NX-1), 2.0/(NY-1), 2.0/(NZ-1)};

    top->SetTopologyType(XDMF_3DCORECTMESH);
    top->GetShapeDesc()->SetShape(nd, cdim);

    geom->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);
    geom->SetOrigin(orig);
    geom->SetDxDyDz(deltas);
#else
    XdmfInt32 nd = 2;
    XdmfInt64 vpdim[] = {nx*NY*NZ, 3};
    XdmfInt64 spdim[] = {nx*NY*NZ, 1};
    XdmfInt64 cdim[] = {(nx-1)*(NY-1)*(NZ-1), 8};

    top->SetTopologyType(XDMF_HEX);
    top->SetNumberOfElements(cdim[0]);

    XdmfArray *top_array = top->GetConnectivity();

#if HEAVY
    sprintf(dsn, "%s.h5:part-%d-connectivity", name, p);
    top_array->SetHeavyDataSetName(dsn);
#endif

    top_array->SetNumberType(XDMF_INT32_TYPE);
    top_array->SetShape(nd, cdim);

    indx = 0;
    for (int i = sx; i < ex; i++)
      for (int j = 0; j < (NY-1); j++)
	for (int k = 0; k < (NZ-1); k++, indx++)
	{
	  int base = (i-sx)*NY*NZ + j*NZ + k;
	  int hex[] = {base, 
	    	   base+1,
		   base+NZ+1,
		   base+NZ,
	  	   base+(NZ*NY), 
	  	   base+(NZ*NY)+1, 
	  	   base+(NZ*NY)+NZ+1, 
	  	   base+(NZ*NY)+NZ};
	  top_array->SetValues(8*indx, hex, 8);
	}

    double orig[] = {-1.0 + (2.0*sx/(NX-1)), -1.0, -1.0};
    double deltas[] = {2.0/(NX-1), 2.0/(NY-1), 2.0/(NZ-1)};

    geom->SetGeometryType(XDMF_GEOMETRY_XYZ);
    XdmfArray *geom_array = geom->GetPoints();

#if HEAVY
    sprintf(dsn, "%s.h5:part-%d-points", name, p);
    geom_array->SetHeavyDataSetName(dsn);
#endif

    geom_array->SetNumberType(XDMF_FLOAT32_TYPE);
    geom_array->SetShape(nd, vpdim);

    indx = 0; float point[3];
    for (int i = sx; i <= ex; i++)
    {
      point[0] = -1.0 + (2.0*i)/(NX-1);
      for (int j = 0; j < NY; j++)
      {
	point[1] = -1.0 + (2.0*j)/(NY-1);
	for (int k = 0; k < NZ; k++, indx++)
	{
	  point[2] = -1.0 + (2.0*k)/(NZ-1);
	  geom_array->SetValues(3*indx, point, 3);
	}
      }
    }

    geom->SetPoints(geom_array);
#endif

    XdmfAttribute *scalar = new XdmfAttribute;
    grid->Insert(scalar);

    scalar->SetName("scalar");
    scalar->SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
    scalar->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);

    XdmfArray *scalar_array = scalar->GetValues();
    scalar_array->SetNumberType(XDMF_FLOAT32_TYPE);

    scalar_array->SetShape(nd, spdim);

    indx = 0;
    for (int i = sx; i <= ex; i++)
    {
      float x = -1.0 + (2.0*i)/(NX-1);
      float xsq = x*x;

      for (int j = 0; j < NY; j++)
      {
	float y = -1.0 + (2.0*j)/(NY-1);
	float ysq = y*y;

	for (int k = 0; k < NZ; k++, indx++)
	{
	  float z = -1.0 + (2.0*k)/(NZ-1);
	  float zsq = z*z;

	  float s = sqrt(xsq + ysq + zsq);
	  scalar_array->SetValue(indx, s);
	}
      }
    }

#if HEAVY
    sprintf(dsn, "%s.h5:part-%d-scalars", name, p);
    scalar_array->SetHeavyDataSetName(dsn);
#endif

    XdmfAttribute *vectors = new XdmfAttribute;
    grid->Insert(vectors);

    vectors->SetName("vectors");
    vectors->SetAttributeType(XDMF_ATTRIBUTE_TYPE_VECTOR);
    vectors->SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);

    XdmfArray *vectors_array = vectors->GetValues();
    vectors_array->SetNumberType(XDMF_FLOAT32_TYPE);

    vectors_array->SetShape(nd, vpdim);

#if HEAVY
    sprintf(dsn, "%s.h5:part-%d-vectors", name, p);
    vectors_array->SetHeavyDataSetName(dsn);
#endif

    XdmfFloat32 vector[] = {0, 0, 0, 1};

    indx = 0;
    for (int i = sx; i <= ex; i++)
    {
      float x = -1.0 + (2.0*i)/(NX-1);
      float xsq = x*x;
      
      vector[1] = dir*x;

      for (int j = 0; j < NY; j++)
      {
	float y = -1.0 + (2.0*j)/(NY-1);
	float ysq = y*y;

	vector[0] = -dir*y;

	for (int k = 0; k < NZ; k++, indx++)
	{
	  float z = -1.0 + (2.0*k)/(NZ-1);
	  float zsq = z*z;

	  float v = sqrt(xsq + ysq + zsq);
	  vectors_array->SetValues(3*indx, vector, 3);
	}
      }
    }
  }

  root.Build();
  sprintf(dsn, "%s.xmf", name);
  ofstream ofs(dsn);
  ofs << root.Serialize() << "\n";
  ofs.close();

  return 1;
}

int
main(int argc, char *argv[])
{
  doit("example",   1);
}
