#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        POWER_LAW
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        2
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  RADIATION                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  Mstar                          0
#define  MU                             1
#define  INJ_TIME                       2
#define  THETA_INJ                      3
#define  D_FLOOR                        4
#define  P_FLOOR                        5
#define  INJ_BUFFER                     6

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MINMOD_LIM
#define  INTERNAL_BOUNDARY              YES
#define  UNIT_LENGTH                    (CONST_au)
#define  UNIT_DENSITY                   3.5e-12
#define  UNIT_VELOCITY                  (sqrt(CONST_G*g_inputParam[Mstar]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))
#define  SHOCK_FLATTENING               MULTID
#define  SHOW_TIME_STEPS                YES
#define  VTK_TIME_INFO                  YES
#define  VTK_VECTOR_DUMP                YES
#define  FAILSALE                       YES
#define  MULTIPLE_LOG_FILES             YES
#define  WARNING_MESSAGES               NO

/* [End] user-defined constants (do not change this line) */
