#include "pluto.h"

void ComputeUserVar(const Data *d, Grid *grid)
/* Created: 14th March 2025, Pawel Janas

PURPOSE: define user-defined output variables (jsu temperature for now)*/
{
    int i, j, k;
    double *x1;

    double ***Te;
    double ***rho, ***prs;
    double ***nH;

    x1 = grid[IDIR].x;

    Te = GetUserVar("Te");
    nH = GetUserVar("nH");

    rho = d->Vc[RHO];
    prs = d->Vc[PRS];

    // const double mu = MeanMolecularWeight(NULL);
    const double mu = g_inputParam[MU];
    DOM_LOOP(k, j, i)
    {
        Te[k][j][i] = prs[k][j][i] / rho[k][j][i] * KELVIN * mu;
        nH[k][j][i] = rho[k][j][i] * UNIT_DENSITY / (mu * CONST_mp);
    }
}


void ChangeDumpVar()
{
  Image *image;

  SetDumpVar("Te", VTK_OUTPUT, YES);
  SetDumpVar("nH", VTK_OUTPUT, YES);
}
