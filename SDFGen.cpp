//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.

#include "makelevelset3.cpp"
#include "config.h"

#include "mex.h"
#include <limits>
#include <cassert>

void mexErrMsgTxt(bool assertion, const char * text)
{
  if(!assertion)
  {
    ::mexErrMsgTxt(text);
  }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  float dx = (float)*mxGetPr(prhs[2]);
  int padding = (int)*mxGetPr(prhs[3]);

  mexErrMsgTxt(nrhs >= 4,"Not enough inputs");
  
  // default
  if(padding < 1) padding = 1;
  //start with a massive inside out bound box.
  Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()), 
    max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());

  const int nv = mxGetM(prhs[0]);
  mexErrMsgTxt(mxGetN(prhs[0]) == 3,"V should be 3D");
  double * V = mxGetPr(prhs[0]);
  const int nf = mxGetM(prhs[1]);
  mexErrMsgTxt(mxGetN(prhs[1]) == 3,"F should be triangles");
  double * F = mxGetPr(prhs[1]);

  std::vector<Vec3f> vertList(nv);
  for(int v = 0;v<nv;v++)
  {
    vertList[v] = Vec3f(V[v+0*nv], V[v+1*nv], V[v+2*nv]);
    update_minmax(vertList[v], min_box, max_box);
  }

  std::vector<Vec3ui> faceList(nf);
  for(int f = 0;f<nf;f++)
  {
    faceList[f] = Vec3ui(F[f+0*nf]-1, F[f+1*nf]-1, F[f+2*nf]-1);
  }

  //Add padding around the box.
  Vec3f unit(1,1,1);
  min_box -= padding*dx*unit;
  max_box += padding*dx*unit;
  Vec3ui sizes = Vec3ui((max_box - min_box)/dx);
  
  Array3f phi_grid;
  make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

  switch(nlhs)
  {
    case 2:
    {
      plhs[1] = mxCreateDoubleMatrix(2,3, mxREAL);
      double * box = mxGetPr(plhs[1]);
      box[0 + 0*2] = min_box[0];
      box[0 + 1*2] = min_box[1];
      box[0 + 2*2] = min_box[2];
      box[1 + 0*2] = max_box[0];
      box[1 + 1*2] = max_box[1];
      box[1 + 2*2] = max_box[2];
      // Fall through
    }
    case 1:
    {
      // Fall through
      std::vector<mwSize> dims = 
        { (mwSize)phi_grid.nj, (mwSize)phi_grid.ni, (mwSize)phi_grid.nk};
      plhs[0] = mxCreateNumericArray(3,&dims[0],mxDOUBLE_CLASS,mxREAL);
      double * phi = mxGetPr(plhs[0]);
      for(unsigned int k = 0; k<dims[2];k++)
      {
        for(unsigned int j = 0; j<dims[1];j++)
        {
          for(unsigned int i = 0; i<dims[0];i++)
          {
            phi[i+dims[0]*(j+dims[1]*k)] = 
              phi_grid.a[j+dims[1]*(i+dims[0]*k)];
          }
        }
      }
    }
    default: break;
  }

}
