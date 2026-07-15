/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author Pawel Janas (pawel.janas@mu.ie)
  \date   July 7th 2026
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "definitions.h"

#define MIN_DENSITY (3.5e-18 / UNIT_DENSITY)

typedef struct
{
  double time_yr;
  double density;
  double pressure;
  double velocity;
} BlastwaveRow;

/* Wrapper for all read in data. */
typedef struct
{
  BlastwaveRow *rows;
  size_t count;
} BlastwaveData;

static BlastwaveData g_blastwave = {.rows = NULL, .count = 0};

BlastwaveData ReadBlastwaveData(const char *filename);
void GetBlastwaveAtTime(double sim_time_yr, double *out_density, double *out_pressure, double *out_velocity);

/* Memory allocation buffers definitions. */
void AllocateShockBuffers(int ni, int nk, int NVARS, int nprocs);

static double *shock_vals = NULL;
static int *shock_found = NULL;
static double *g_shock_vals = NULL;
static int *g_shock_found = NULL;

static int buffers_initialised = 0;
/* --------------------------------------------------------------------- */

/* ********************************************************************* */
void Init(double *v, double x1, double x2, double x3)
/*!
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical}
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  g_minCoolingTemp = 12.0;                 // limit cooling before blastwave injection
  g_maxCoolingRate = 0.9;                  // default is 0.1
  g_smallPressure = g_inputParam[P_FLOOR]; // default is 1e-12
}

/* ********************************************************************* */
void InitDomain(Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  /* Load in 1D blastwave data. */
  int static do_once = 1;
  if (do_once)
  {
    const char *filename = "1D_blastwave_Close_extend.txt";
    g_blastwave = ReadBlastwaveData(filename);
    if (g_blastwave.count == 0)
    {
      fprintf(stderr, "Error: No data read from %s\n", filename);
      exit(EXIT_FAILURE);
    }
    printf("Blastwave data loaded: %zu rows.\n", g_blastwave.count);
    do_once = 0;
  }

  /* LOAD STABLE DISK */
  int i, j, k, id, size[3];
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  long int offset;

  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 0, CENTER);
  if (id < 0)
  {
    fprintf(stderr, "Error opening input data file.\n");
    exit(EXIT_FAILURE);
  }

  InputDataGridSize(id, size);
  offset = (long)size[0] * (long)size[1] * (long)size[2];

  /* Density */
  TOT_LOOP(k, j, i)
  {
    d->Vc[RHO][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  /* Velocity */
  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 1 * offset, CENTER); // v_r (VX1)
  TOT_LOOP(k, j, i)
  {
    d->Vc[VX1][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 2 * offset, CENTER); // v_phi (VX2)
  TOT_LOOP(k, j, i)
  {
    d->Vc[VX2][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 3 * offset, CENTER); // v_theta (VX3)
  TOT_LOOP(k, j, i)
  {
    d->Vc[VX3][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  /* Pressure */
  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 4 * offset, CENTER);
  TOT_LOOP(k, j, i)
  {
    d->Vc[PRS][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  /* Tracer (TRC)*/
  id = InputDataOpen("INPUT_DATA/highmu/data.1000.dbl", "INPUT_DATA/highmu/grid.out", " ", 5 * offset, CENTER);
  TOT_LOOP(k, j, i)
  {
    d->Vc[TRC][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
  }
  InputDataClose(id);

  /* Set disk tracer (TRC+1 = 0)*/
  TOT_LOOP(k, j, i)
  {
    if (d->Vc[TRC][k][j][i] < 1.005)
      d->Vc[TRC + 1][k][j][i] = 0.0;
  }
}

/* ********************************************************************* */
void Analysis(const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField(double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
  B0[0] = 0.0;
  B0[1] = 0.0;
  B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary(const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END,
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  double *x1, *x2, *x3, R, r, Omega_k, Omega;
  double h = 0.05; // disk aspect ratio
  int do_once = 1;
  static int forced_small_dt_done = 0;

  double T_amb = 1e4;   // ambient temperature (Kelvin)
  double mu_amb = 0.71; // ambient ean molecular weight
  // double mu = MeanMolecularWeight(d->Vc);

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  // Define injection parameters
  double sim_time = g_time;
  double sim_time_yr = sim_time * (UNIT_LENGTH / UNIT_VELOCITY) / 31557600; // convert simulation time to years (Julian astronomical year).
  double theta_inj = g_inputParam[THETA_INJ] * CONST_PI / 180;              // convert input direction parameter into radians
  double r_threshold = 140.0;
  double r_threshold_perp = 400.0; // threshold for flagging cells perpendicular to injection slab
  double slab_thickness = 20.0;

  /* Enforce a small dt at blastwave injection time to prevent abort. */
  if (sim_time_yr >= g_inputParam[INJ_TIME] && !forced_small_dt_done)
  {
    g_dt = 1e-8;
    forced_small_dt_done = 1;
  }

/*----------------------------------------------------------------------------------------*/  
/*-----------------------BLASTWAVE / WIND INJECTION---------------------------------------*/
/*----------------------------------------------------------------------------------------*/  
  /* For constant flow (wind). */
  double UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  double rho_blast = 5.74e-21 / UNIT_DENSITY;
  double vel_blast = 2.87 * (1e8 / UNIT_VELOCITY);
  double P_blast;

  /* For dynamic blastwave. Interpolate the read in file at Init(). */
  // double UNIT_PRESSURE = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  // double rho_blast, vel_blast, P_blast;
  // double inj_time_yr = sim_time_yr - g_inputParam[INJ_TIME]; // make sure blastwave starts at INJ_TIME.
  // inj_time_yr += g_inputParam[INJ_BUFFER];                   // to induce faster interaction with blastwave
  // if (inj_time_yr >= 0)
  // {
  //   GetBlastwaveAtTime(inj_time_yr, &rho_blast, &P_blast, &vel_blast);
  // }
  rho_blast *= (3.5e-18 / UNIT_DENSITY); // Input data is in units of 3.5e-18 g cm-3. Converting to 3.5e-12 units.
  P_blast *= (3.5e-2 / UNIT_PRESSURE);   // Input data is in units of 3.5e-2 dyne cm-2. Converting to code units.
  vel_blast *= (1e8 / UNIT_VELOCITY);    // Input data is in units of 1e8 cm/s. Converting to code units.

  // #if DIMENSIONS == 3
  if (side == 0)
  {
    if (sim_time_yr >= g_inputParam[INJ_TIME])
    {
      g_minCoolingTemp = 12.0; // allow cooling to lower temperatures after blastwave injection

      DOM_LOOP(k, j, i)
      {
        r = x1[i];
        double x = r * cos(x2[j]);
        double y = r * sin(x2[j]);
        double slab_coord = x * cos(theta_inj) + y * sin(theta_inj);

        /* Flag cells perpendicular to injection slab. */
        double angle_perp1 = theta_inj + CONST_PI / 2;
        double angle_perp2 = theta_inj - CONST_PI / 2;
        double angle_perp2_offset = angle_perp2 - (CONST_PI / 32); // offset to allow more cells in the perpendicular direction near the injection region
        double angle_perp1_offset = angle_perp1 + (CONST_PI / 32); // offset to allow more cells in the perpendicular direction near the injection region
        double slab_coord_perp1 = x * cos(angle_perp1) + y * sin(angle_perp1);
        double slab_coord_perp2 = x * cos(angle_perp2) + y * sin(angle_perp2);
        double slab_coord_perp2_offset = x * cos(angle_perp2_offset) + y * sin(angle_perp2_offset);
        double slab_coord_perp1_offset = x * cos(angle_perp1_offset) + y * sin(angle_perp1_offset);

        if ((slab_coord_perp1 > r_threshold_perp && slab_coord_perp1_offset > r_threshold_perp + 26) || (slab_coord_perp2 > r_threshold_perp && slab_coord_perp2_offset > r_threshold_perp + 26))
        {
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY; // don't simulate
          continue;
        }

        if (slab_coord > r_threshold) // && (slab_coord_perp2 < r_threshold || slab_coord_perp1 < r_threshold))
        {
          d->Vc[RHO][k][j][i] = rho_blast;
          d->Vc[VX1][k][j][i] = -vel_blast * cos(theta_inj - x2[j]);
          d->Vc[VX2][k][j][i] = -vel_blast * sin(theta_inj - x2[j]);
          d->Vc[PRS][k][j][i] = P_blast;
          d->Vc[TRC + 1][k][j][i] = 1.0; // tracer for blastwave/wind material
          // d->Vc[PRS][k][j][i] = T_blast * d->Vc[RHO][k][j][i] / (KELVIN * mu_amb);
        }
/*----------------------------------------------------------------------------------------*/  
/*----------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------*/  

        double dth = fabs(x2[j] - theta_inj); // angular distance from injection direction
        double dth_pole = fabs(x2[j] - CONST_PI);
        double dth_pole_0 = fabs(x2[j]);

        // Handle wraparound for θ = 0 (cells at θ ≈ 2π)
        if (dth_pole_0 > CONST_PI)
        {
          dth_pole_0 = 2.0 * CONST_PI - dth_pole_0;
        }
        if (dth_pole > CONST_PI)
        {
          dth_pole = 2.0 * CONST_PI - dth_pole;
        }
        if (dth > CONST_PI)
        {
          dth = 2.0 * CONST_PI - dth;
        }

        // Near shock front: 0.91 rad ~ 52 degrees (based on visual analysis).
        if (dth < 0.91)
        {
          d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], P_blast);
          d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], rho_blast);
        }
        // Wider arc to try get plane-parallel blastwave acting on disk - prevents weaker regions near turned off cells.
        // Update 15/04/26: Impose blastwave velocities (no MAX call). - this didnt work
        // Update 16/04/26: Enforce the blastwave prs and rho values in the wider arc. - this also didn't work.
        // Update 16/04/26 (2): Move pressure and density towards blastwave values, limiting to 10% change per timestep. - this didn't work.
        // update 16/04/26 (3): Enforce blastwave values in wider arc, but only if they are above the floor.
        // Update 17/04/26: This is working well, but at early times, the blastwave shocked cells are forced to the blastwave values.
        // Need to implement a changing radius threshold depending on blastwave strength.
        else if (r > r_threshold_perp && dth < 1.32) // 1.32 rad ~ 75 degrees
        {
          // d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], P_blast);
          // d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], rho_blast);
          double p_state = d->Vc[PRS][k][j][i];
          double rho_state = d->Vc[RHO][k][j][i];

          double p_target = P_blast;
          double rho_target = rho_blast;

          double decay_alpha = 0.1; // smaller = max holds longer
          if (p_target > p_state)
          {
            p_state = p_target;
          }
          else
          {
            p_state = p_state * (1.0 - decay_alpha) + p_target * decay_alpha;
          }
          if (rho_target > rho_state)
          {
            rho_state = rho_target;
          }
          else
          {
            rho_state = rho_state * (1.0 - decay_alpha) + rho_target * decay_alpha;
          }
          d->Vc[PRS][k][j][i] = MAX(p_state, g_inputParam[P_FLOOR]);
          d->Vc[RHO][k][j][i] = MAX(rho_state, g_inputParam[D_FLOOR]);
        }
        else
        {
          d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], g_inputParam[P_FLOOR]);
          d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], g_inputParam[D_FLOOR]);
        }
      }
    }
    else
    {
      DOM_LOOP(k, j, i)
      {
        r = x1[i];
        d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], g_inputParam[D_FLOOR]);
        d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], g_inputParam[P_FLOOR]);
      }
    }
  }
  // #endif

  double v[256];
  if (side == X1_BEG)
  {
    X1_BEG_LOOP(k, j, i)
    {
      // Init(v, x1[i], x2[j], x3[k]);
      //
      // d->Vc[nv][k][j][i] = v[nv]; // copy values from Init()
      NVAR_LOOP(nv)
      d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2 * IBEG - i - 1];

      /* Validate cell values and enforce floors */
      d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], g_inputParam[D_FLOOR]);
      d->Vc[PRS][k][j][i] = MAX(d->Vc[PRS][k][j][i], g_inputParam[P_FLOOR]);

      if (!isfinite(d->Vc[VX1][k][j][i]))
        d->Vc[VX1][k][j][i] = 0.0;
      if (!isfinite(d->Vc[VX2][k][j][i]))
        d->Vc[VX2][k][j][i] = 0.0;
      if (!isfinite(d->Vc[VX3][k][j][i]))
        d->Vc[VX3][k][j][i] = 0.0;

      /* Diode BC, no inflow */
      if (!isfinite(d->Vc[VX1][k][j][i]) || d->Vc[VX1][k][j][i] > 0.0)
        d->Vc[VX1][k][j][i] = 0.0;

#if GEOMETRY == SPHERICAL
      r = x1[i];
#if DIMENSIONS == 2
      R = x1[i] * sin(x2[j]);
#elif DIMENSIONS == 3
      R = x1[i] * sin(x2[j]);
#endif
#endif
      /* softening as in Init(): use R_soft to compute Omega_k */
      const double R_floor = 1e-6;
      double R_cyl = fabs(R);
      if (R_cyl < R_floor)
        R_cyl = R_floor;
      const double R_soft_eps = 1e-9;
      double R_soft = sqrt(R_cyl * R_cyl + R_soft_eps * R_soft_eps);

      /* At poles or very small R, set velocities to zero */
      if (R_cyl < 1e-5 || sin(x2[j]) < 1e-4)
      {
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[VX3][k][j][i] = 0.0;
      }
      else
      {
        Omega_k = 2.0 * CONST_PI / (R_soft * sqrt(R_soft));
        double omega_arg = R_cyl / r - 5.0 / 2.0 * h * h;
        if (omega_arg < 0.0)
          omega_arg = 0.0;
        Omega = Omega_k * sqrt(omega_arg);
        if (!isfinite(Omega) || !isfinite(R))
          d->Vc[VX3][k][j][i] = 0.0;
        else
          d->Vc[VX3][k][j][i] = R_cyl * Omega;
      }
    }
  }

  if (side == X1_END)
  {
    X1_END_LOOP(k, j, i)
    {
      NVAR_LOOP(nv)
      d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      if (d->Vc[VX1][k][j][i] < 0.0)
        d->Vc[VX1][k][j][i] = 0.0; // zero-gradient outflow BC (prevent inflow)
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = -4.0 * CONST_PI * CONST_PI / x1 / x1;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  double R, r, z, th, phi;
  r = x1;
  // R = x1 * sin(x2); // cylindrical radius

  phi = -4.0 * CONST_PI * CONST_PI / r;

  return phi;
}
#endif

BlastwaveData ReadBlastwaveData(const char *filename)
{
  FILE *file = fopen(filename, "r");
  if (!file)
  {
    fprintf(stderr, "Error opening file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  char line[1024];
  size_t capacity = 200; // initial guess for capacity
  size_t count = 0;

  BlastwaveRow *rows = malloc(capacity * sizeof(BlastwaveRow));
  if (!rows)
  {
    fprintf(stderr, "Memory allocation failed\n");
    fclose(file);
    exit(EXIT_FAILURE);
  }

  if (fgets(line, sizeof(line), file) == NULL)
  {
    fprintf(stderr, "Error reading header from file %s\n", filename);
    fclose(file);
    free(rows);
    exit(EXIT_FAILURE);
  }

  while (fgets(line, sizeof(line), file) != NULL)
  {
    if (count == capacity)
    {
      capacity *= 2; // double the capacity
      BlastwaveRow *tmp = realloc(rows, capacity * sizeof(BlastwaveRow));
      if (!tmp)
      {
        fprintf(stderr, "Memory reallocation failed\n");
        fclose(file);
        free(rows);
        exit(EXIT_FAILURE);
      }
      rows = tmp;
    }

    int n = sscanf(line, "%lf\t%lf\t%lf\t%lf",
                   &rows[count].time_yr,
                   &rows[count].density,
                   &rows[count].pressure,
                   &rows[count].velocity);

    if (n == 4)
    {
      count++;
    }
    else
    {
      fprintf(stderr, "Skipping malformed line: %s\n", line);
    }
  }

  fclose(file);

  /* Shrink arrays to exact size. */
  if (count < capacity)
  {
    BlastwaveRow *tmp = realloc(rows, count * sizeof(BlastwaveRow));
    if (tmp)
    {
      rows = tmp;
    }
  }

  return (BlastwaveData){.rows = rows, .count = count};
}

void GetBlastwaveAtTime(double sim_time_yr, double *out_density, double *out_pressure, double *out_velocity)
/*!
 * Given a simulation time sim_time_yr, look up (and linearly interpolate)
 * density, pressure, velocity from the global table g_blastwave.
 *
 *   sim_time_yr : current simulation time [yr] -- IMPORTANT!!!
 *   out_density  : pointer where to store the interpolated density
 *   out_pressure : pointer where to store the interpolated pressure
 *   out_velocity : pointer where to store the interpolated velocity
 *
 * If sim_time_yr is:
 *   - below the first row’s time, returns the first row’s values.
 *   - above the last row’s time, returns the last row’s values.
 *   - in between, returns the linearly‐interpolated values.
 */
{
  size_t n = g_blastwave.count;
  if (n == 0)
  {
    fprintf(stderr, "Error: No blastwave data available.\n");
    exit(EXIT_FAILURE);
  }

  BlastwaveRow *R = g_blastwave.rows;

  /* If before first item, clamp to row 0. */
  if (sim_time_yr <= R[0].time_yr)
  {
    *out_density = R[0].density;
    *out_pressure = R[0].pressure;
    *out_velocity = R[0].velocity;
    return;
  }

  /* If after last item, clamp to last row. */
  if (sim_time_yr >= R[n - 1].time_yr)
  {
    *out_density = R[n - 1].density;
    *out_pressure = R[n - 1].pressure;
    *out_velocity = R[n - 1].velocity;
    return;
  }

  /* Else, find i such that R[i].time <= sim_time_yr <= R[i+1].time. */
  size_t i_lo = 0, i_hi = n - 1;
  /* linear search through rows for desired sim_time. */
  for (size_t i = 0; i + 1 < n; i++)
  {
    if (R[i].time_yr <= sim_time_yr && sim_time_yr <= R[i + 1].time_yr)
    {
      i_lo = i;
      i_hi = i + 1;
      break;
    }
  }

  double t_lo = R[i_lo].time_yr;
  double t_hi = R[i_hi].time_yr;

  /* Avoid divide by zero. */
  double alpha;
  if (t_hi == t_lo)
  {
    alpha = 0.0;
  }
  else
  {
    alpha = (sim_time_yr - t_lo) / (t_hi - t_lo);
  }

  /* Linearly interpolate each quantity. */
  *out_density = R[i_lo].density + alpha * (R[i_hi].density - R[i_lo].density);
  *out_pressure = R[i_lo].pressure + alpha * (R[i_hi].pressure - R[i_lo].pressure);
  *out_velocity = R[i_lo].velocity + alpha * (R[i_hi].velocity - R[i_lo].velocity);
}

void AllocateShockBuffers(int ni, int nk, int NVARS, int nprocs)
{
  if (buffers_initialised)
    return;

  int local_cells = ni * nk;

  shock_vals = (double *)calloc((size_t)local_cells * (size_t)NVARS, sizeof(double));
  shock_found = (int *)calloc((size_t)local_cells, sizeof(int));
  g_shock_vals = (double *)calloc((size_t)local_cells * (size_t)NVARS, sizeof(double));
  g_shock_found = (int *)calloc((size_t)local_cells, sizeof(int));
  // send_idx = malloc(local_cells * sizeof(int));
  // recv_idx = malloc(global_cells * sizeof(int));
  // recv_counts = malloc(nprocs * sizeof(int));
  // recv_displs = malloc(nprocs * sizeof(int));
  // send_vals = malloc(local_cells * NVARS * sizeof(double));
  // recv_vals = malloc(global_cells * NVARS * sizeof(double));

  // recv_counts_vals = malloc(nprocs * sizeof(int));
  // recv_displs_vals = malloc(nprocs * sizeof(int));

  if (!shock_vals || !shock_found || !g_shock_vals || !g_shock_found)
  {
    fprintf(stderr, "Error: Shock buffer allocation failed\n");
    exit(EXIT_FAILURE);
  }

  buffers_initialised = 1;
}

// void FreeShockBuffers(void)
// {
//   if (!buffers_initialised)
//     return;

//   free(shock_vals);
//   free(shock_found);
//   free(g_shock_vals);
//   free(g_shock_found);
//   free(send_idx);
//   free(recv_idx);
//   free(recv_counts);
//   free(recv_displs);
//   free(send_vals);
//   free(recv_vals);
//   free(recv_counts_vals);
//   free(recv_displs_vals);

//   shock_vals = NULL;
//   shock_found = NULL;
//   g_shock_vals = NULL;
//   g_shock_found = NULL;
//   send_idx = NULL;
//   recv_idx = NULL;
//   recv_counts = NULL;
//   recv_displs = NULL;
//   send_vals = NULL;
//   recv_vals = NULL;
//   recv_counts_vals = NULL;
//   recv_displs_vals = NULL;

//   buffers_initialised = 0;
// }