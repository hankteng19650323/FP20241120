#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7049607932794873382) {
   out_7049607932794873382[0] = delta_x[0] + nom_x[0];
   out_7049607932794873382[1] = delta_x[1] + nom_x[1];
   out_7049607932794873382[2] = delta_x[2] + nom_x[2];
   out_7049607932794873382[3] = delta_x[3] + nom_x[3];
   out_7049607932794873382[4] = delta_x[4] + nom_x[4];
   out_7049607932794873382[5] = delta_x[5] + nom_x[5];
   out_7049607932794873382[6] = delta_x[6] + nom_x[6];
   out_7049607932794873382[7] = delta_x[7] + nom_x[7];
   out_7049607932794873382[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_714116200904137688) {
   out_714116200904137688[0] = -nom_x[0] + true_x[0];
   out_714116200904137688[1] = -nom_x[1] + true_x[1];
   out_714116200904137688[2] = -nom_x[2] + true_x[2];
   out_714116200904137688[3] = -nom_x[3] + true_x[3];
   out_714116200904137688[4] = -nom_x[4] + true_x[4];
   out_714116200904137688[5] = -nom_x[5] + true_x[5];
   out_714116200904137688[6] = -nom_x[6] + true_x[6];
   out_714116200904137688[7] = -nom_x[7] + true_x[7];
   out_714116200904137688[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8382138516203153115) {
   out_8382138516203153115[0] = 1.0;
   out_8382138516203153115[1] = 0;
   out_8382138516203153115[2] = 0;
   out_8382138516203153115[3] = 0;
   out_8382138516203153115[4] = 0;
   out_8382138516203153115[5] = 0;
   out_8382138516203153115[6] = 0;
   out_8382138516203153115[7] = 0;
   out_8382138516203153115[8] = 0;
   out_8382138516203153115[9] = 0;
   out_8382138516203153115[10] = 1.0;
   out_8382138516203153115[11] = 0;
   out_8382138516203153115[12] = 0;
   out_8382138516203153115[13] = 0;
   out_8382138516203153115[14] = 0;
   out_8382138516203153115[15] = 0;
   out_8382138516203153115[16] = 0;
   out_8382138516203153115[17] = 0;
   out_8382138516203153115[18] = 0;
   out_8382138516203153115[19] = 0;
   out_8382138516203153115[20] = 1.0;
   out_8382138516203153115[21] = 0;
   out_8382138516203153115[22] = 0;
   out_8382138516203153115[23] = 0;
   out_8382138516203153115[24] = 0;
   out_8382138516203153115[25] = 0;
   out_8382138516203153115[26] = 0;
   out_8382138516203153115[27] = 0;
   out_8382138516203153115[28] = 0;
   out_8382138516203153115[29] = 0;
   out_8382138516203153115[30] = 1.0;
   out_8382138516203153115[31] = 0;
   out_8382138516203153115[32] = 0;
   out_8382138516203153115[33] = 0;
   out_8382138516203153115[34] = 0;
   out_8382138516203153115[35] = 0;
   out_8382138516203153115[36] = 0;
   out_8382138516203153115[37] = 0;
   out_8382138516203153115[38] = 0;
   out_8382138516203153115[39] = 0;
   out_8382138516203153115[40] = 1.0;
   out_8382138516203153115[41] = 0;
   out_8382138516203153115[42] = 0;
   out_8382138516203153115[43] = 0;
   out_8382138516203153115[44] = 0;
   out_8382138516203153115[45] = 0;
   out_8382138516203153115[46] = 0;
   out_8382138516203153115[47] = 0;
   out_8382138516203153115[48] = 0;
   out_8382138516203153115[49] = 0;
   out_8382138516203153115[50] = 1.0;
   out_8382138516203153115[51] = 0;
   out_8382138516203153115[52] = 0;
   out_8382138516203153115[53] = 0;
   out_8382138516203153115[54] = 0;
   out_8382138516203153115[55] = 0;
   out_8382138516203153115[56] = 0;
   out_8382138516203153115[57] = 0;
   out_8382138516203153115[58] = 0;
   out_8382138516203153115[59] = 0;
   out_8382138516203153115[60] = 1.0;
   out_8382138516203153115[61] = 0;
   out_8382138516203153115[62] = 0;
   out_8382138516203153115[63] = 0;
   out_8382138516203153115[64] = 0;
   out_8382138516203153115[65] = 0;
   out_8382138516203153115[66] = 0;
   out_8382138516203153115[67] = 0;
   out_8382138516203153115[68] = 0;
   out_8382138516203153115[69] = 0;
   out_8382138516203153115[70] = 1.0;
   out_8382138516203153115[71] = 0;
   out_8382138516203153115[72] = 0;
   out_8382138516203153115[73] = 0;
   out_8382138516203153115[74] = 0;
   out_8382138516203153115[75] = 0;
   out_8382138516203153115[76] = 0;
   out_8382138516203153115[77] = 0;
   out_8382138516203153115[78] = 0;
   out_8382138516203153115[79] = 0;
   out_8382138516203153115[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_1977712703250497707) {
   out_1977712703250497707[0] = state[0];
   out_1977712703250497707[1] = state[1];
   out_1977712703250497707[2] = state[2];
   out_1977712703250497707[3] = state[3];
   out_1977712703250497707[4] = state[4];
   out_1977712703250497707[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1977712703250497707[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1977712703250497707[7] = state[7];
   out_1977712703250497707[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6944128836838434191) {
   out_6944128836838434191[0] = 1;
   out_6944128836838434191[1] = 0;
   out_6944128836838434191[2] = 0;
   out_6944128836838434191[3] = 0;
   out_6944128836838434191[4] = 0;
   out_6944128836838434191[5] = 0;
   out_6944128836838434191[6] = 0;
   out_6944128836838434191[7] = 0;
   out_6944128836838434191[8] = 0;
   out_6944128836838434191[9] = 0;
   out_6944128836838434191[10] = 1;
   out_6944128836838434191[11] = 0;
   out_6944128836838434191[12] = 0;
   out_6944128836838434191[13] = 0;
   out_6944128836838434191[14] = 0;
   out_6944128836838434191[15] = 0;
   out_6944128836838434191[16] = 0;
   out_6944128836838434191[17] = 0;
   out_6944128836838434191[18] = 0;
   out_6944128836838434191[19] = 0;
   out_6944128836838434191[20] = 1;
   out_6944128836838434191[21] = 0;
   out_6944128836838434191[22] = 0;
   out_6944128836838434191[23] = 0;
   out_6944128836838434191[24] = 0;
   out_6944128836838434191[25] = 0;
   out_6944128836838434191[26] = 0;
   out_6944128836838434191[27] = 0;
   out_6944128836838434191[28] = 0;
   out_6944128836838434191[29] = 0;
   out_6944128836838434191[30] = 1;
   out_6944128836838434191[31] = 0;
   out_6944128836838434191[32] = 0;
   out_6944128836838434191[33] = 0;
   out_6944128836838434191[34] = 0;
   out_6944128836838434191[35] = 0;
   out_6944128836838434191[36] = 0;
   out_6944128836838434191[37] = 0;
   out_6944128836838434191[38] = 0;
   out_6944128836838434191[39] = 0;
   out_6944128836838434191[40] = 1;
   out_6944128836838434191[41] = 0;
   out_6944128836838434191[42] = 0;
   out_6944128836838434191[43] = 0;
   out_6944128836838434191[44] = 0;
   out_6944128836838434191[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6944128836838434191[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6944128836838434191[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6944128836838434191[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6944128836838434191[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6944128836838434191[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6944128836838434191[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6944128836838434191[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6944128836838434191[53] = -9.8000000000000007*dt;
   out_6944128836838434191[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6944128836838434191[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6944128836838434191[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6944128836838434191[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6944128836838434191[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6944128836838434191[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6944128836838434191[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6944128836838434191[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6944128836838434191[62] = 0;
   out_6944128836838434191[63] = 0;
   out_6944128836838434191[64] = 0;
   out_6944128836838434191[65] = 0;
   out_6944128836838434191[66] = 0;
   out_6944128836838434191[67] = 0;
   out_6944128836838434191[68] = 0;
   out_6944128836838434191[69] = 0;
   out_6944128836838434191[70] = 1;
   out_6944128836838434191[71] = 0;
   out_6944128836838434191[72] = 0;
   out_6944128836838434191[73] = 0;
   out_6944128836838434191[74] = 0;
   out_6944128836838434191[75] = 0;
   out_6944128836838434191[76] = 0;
   out_6944128836838434191[77] = 0;
   out_6944128836838434191[78] = 0;
   out_6944128836838434191[79] = 0;
   out_6944128836838434191[80] = 1;
}
void h_25(double *state, double *unused, double *out_3050096150093046930) {
   out_3050096150093046930[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4007907610137106257) {
   out_4007907610137106257[0] = 0;
   out_4007907610137106257[1] = 0;
   out_4007907610137106257[2] = 0;
   out_4007907610137106257[3] = 0;
   out_4007907610137106257[4] = 0;
   out_4007907610137106257[5] = 0;
   out_4007907610137106257[6] = 1;
   out_4007907610137106257[7] = 0;
   out_4007907610137106257[8] = 0;
}
void h_24(double *state, double *unused, double *out_310778294690262708) {
   out_310778294690262708[0] = state[4];
   out_310778294690262708[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8754416937103615291) {
   out_8754416937103615291[0] = 0;
   out_8754416937103615291[1] = 0;
   out_8754416937103615291[2] = 0;
   out_8754416937103615291[3] = 0;
   out_8754416937103615291[4] = 1;
   out_8754416937103615291[5] = 0;
   out_8754416937103615291[6] = 0;
   out_8754416937103615291[7] = 0;
   out_8754416937103615291[8] = 0;
   out_8754416937103615291[9] = 0;
   out_8754416937103615291[10] = 0;
   out_8754416937103615291[11] = 0;
   out_8754416937103615291[12] = 0;
   out_8754416937103615291[13] = 0;
   out_8754416937103615291[14] = 1;
   out_8754416937103615291[15] = 0;
   out_8754416937103615291[16] = 0;
   out_8754416937103615291[17] = 0;
}
void h_30(double *state, double *unused, double *out_4570068002176880642) {
   out_4570068002176880642[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1489574651629857630) {
   out_1489574651629857630[0] = 0;
   out_1489574651629857630[1] = 0;
   out_1489574651629857630[2] = 0;
   out_1489574651629857630[3] = 0;
   out_1489574651629857630[4] = 1;
   out_1489574651629857630[5] = 0;
   out_1489574651629857630[6] = 0;
   out_1489574651629857630[7] = 0;
   out_1489574651629857630[8] = 0;
}
void h_26(double *state, double *unused, double *out_5929695533186250470) {
   out_5929695533186250470[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7749410929011162481) {
   out_7749410929011162481[0] = 0;
   out_7749410929011162481[1] = 0;
   out_7749410929011162481[2] = 0;
   out_7749410929011162481[3] = 0;
   out_7749410929011162481[4] = 0;
   out_7749410929011162481[5] = 0;
   out_7749410929011162481[6] = 0;
   out_7749410929011162481[7] = 1;
   out_7749410929011162481[8] = 0;
}
void h_27(double *state, double *unused, double *out_9152615585826241320) {
   out_9152615585826241320[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3664337963430282541) {
   out_3664337963430282541[0] = 0;
   out_3664337963430282541[1] = 0;
   out_3664337963430282541[2] = 0;
   out_3664337963430282541[3] = 1;
   out_3664337963430282541[4] = 0;
   out_3664337963430282541[5] = 0;
   out_3664337963430282541[6] = 0;
   out_3664337963430282541[7] = 0;
   out_3664337963430282541[8] = 0;
}
void h_29(double *state, double *unused, double *out_288504160933220134) {
   out_288504160933220134[0] = state[1];
}
void H_29(double *state, double *unused, double *out_979343307315465446) {
   out_979343307315465446[0] = 0;
   out_979343307315465446[1] = 1;
   out_979343307315465446[2] = 0;
   out_979343307315465446[3] = 0;
   out_979343307315465446[4] = 0;
   out_979343307315465446[5] = 0;
   out_979343307315465446[6] = 0;
   out_979343307315465446[7] = 0;
   out_979343307315465446[8] = 0;
}
void h_28(double *state, double *unused, double *out_3651398616802353461) {
   out_3651398616802353461[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6061742324384996020) {
   out_6061742324384996020[0] = 1;
   out_6061742324384996020[1] = 0;
   out_6061742324384996020[2] = 0;
   out_6061742324384996020[3] = 0;
   out_6061742324384996020[4] = 0;
   out_6061742324384996020[5] = 0;
   out_6061742324384996020[6] = 0;
   out_6061742324384996020[7] = 0;
   out_6061742324384996020[8] = 0;
}
void h_31(double *state, double *unused, double *out_6258930638932653504) {
   out_6258930638932653504[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8375619031244513957) {
   out_8375619031244513957[0] = 0;
   out_8375619031244513957[1] = 0;
   out_8375619031244513957[2] = 0;
   out_8375619031244513957[3] = 0;
   out_8375619031244513957[4] = 0;
   out_8375619031244513957[5] = 0;
   out_8375619031244513957[6] = 0;
   out_8375619031244513957[7] = 0;
   out_8375619031244513957[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_7049607932794873382) {
  err_fun(nom_x, delta_x, out_7049607932794873382);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_714116200904137688) {
  inv_err_fun(nom_x, true_x, out_714116200904137688);
}
void car_H_mod_fun(double *state, double *out_8382138516203153115) {
  H_mod_fun(state, out_8382138516203153115);
}
void car_f_fun(double *state, double dt, double *out_1977712703250497707) {
  f_fun(state,  dt, out_1977712703250497707);
}
void car_F_fun(double *state, double dt, double *out_6944128836838434191) {
  F_fun(state,  dt, out_6944128836838434191);
}
void car_h_25(double *state, double *unused, double *out_3050096150093046930) {
  h_25(state, unused, out_3050096150093046930);
}
void car_H_25(double *state, double *unused, double *out_4007907610137106257) {
  H_25(state, unused, out_4007907610137106257);
}
void car_h_24(double *state, double *unused, double *out_310778294690262708) {
  h_24(state, unused, out_310778294690262708);
}
void car_H_24(double *state, double *unused, double *out_8754416937103615291) {
  H_24(state, unused, out_8754416937103615291);
}
void car_h_30(double *state, double *unused, double *out_4570068002176880642) {
  h_30(state, unused, out_4570068002176880642);
}
void car_H_30(double *state, double *unused, double *out_1489574651629857630) {
  H_30(state, unused, out_1489574651629857630);
}
void car_h_26(double *state, double *unused, double *out_5929695533186250470) {
  h_26(state, unused, out_5929695533186250470);
}
void car_H_26(double *state, double *unused, double *out_7749410929011162481) {
  H_26(state, unused, out_7749410929011162481);
}
void car_h_27(double *state, double *unused, double *out_9152615585826241320) {
  h_27(state, unused, out_9152615585826241320);
}
void car_H_27(double *state, double *unused, double *out_3664337963430282541) {
  H_27(state, unused, out_3664337963430282541);
}
void car_h_29(double *state, double *unused, double *out_288504160933220134) {
  h_29(state, unused, out_288504160933220134);
}
void car_H_29(double *state, double *unused, double *out_979343307315465446) {
  H_29(state, unused, out_979343307315465446);
}
void car_h_28(double *state, double *unused, double *out_3651398616802353461) {
  h_28(state, unused, out_3651398616802353461);
}
void car_H_28(double *state, double *unused, double *out_6061742324384996020) {
  H_28(state, unused, out_6061742324384996020);
}
void car_h_31(double *state, double *unused, double *out_6258930638932653504) {
  h_31(state, unused, out_6258930638932653504);
}
void car_H_31(double *state, double *unused, double *out_8375619031244513957) {
  H_31(state, unused, out_8375619031244513957);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
