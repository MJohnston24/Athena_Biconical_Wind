//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

void JetInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void cgbc_i(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void zero_mom_i(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

// ambient density and pressure and E0 sent from the input file
static Real gam;
static Real pa;
// static Real prat;
static Real rhoa;
static Real E0;
static Real Mach;
static Real vr;
static Real rs;
static Real theta1;
static Real theta2;
static Real t_stop;
static Real rmin;

Real shockProperties(MeshBlock *pmb,int iout);

void Mesh::InitUserMeshData(ParameterInput *pin) {

  if (pin->GetString("mesh","ix1_bc") == "user") {
    // enroll boundary value function pointers
    std::cout << "Enrolling continuous energy injection" << std::endl;
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, JetInnerX1);
  }


  // AllocateUserHistoryOutput(6);
  // EnrollUserHistoryOutput(0, shockProperties, "Esh");
  // EnrollUserHistoryOutput(1, shockProperties, "Rsh");
  // EnrollUserHistoryOutput(2, shockProperties, "Thetash");
  // EnrollUserHistoryOutput(3, shockProperties, "Psh");
  // EnrollUserHistoryOutput(4, shockProperties, "Rhosh");
  // EnrollUserHistoryOutput(5, shockProperties, "Ush");
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  pa      = pin->GetOrAddReal("problem", "pamb", 1.0);
  // prat    = pin->GetOrAddReal("problem", "prat", 1.0);
  rhoa    = pin->GetOrAddReal("problem", "rhoamb", 1.0);
  E0      = pin->GetReal("problem", "E0");
  Mach    = pin->GetReal("problem", "Mach");
  vr      = pin->GetReal("problem", "vr");
  rs      = pin->GetReal("problem", "rs");
  theta1  = pin->GetReal("problem", "theta1");
  theta2  = pin->GetReal("problem", "theta2");
  t_stop  = pin->GetReal("problem", "tstop");
  rmin    = pin->GetReal("mesh", "x1min");

  gam = peos->GetGamma();
  Real gm1 = gam - 1.0;
  // Real rsh = pcoord->x1f(is+1);  // this gives xi = 1.33
  Real rsh = pcoord->x1v(is); // this gives xi = 1.23
  // volume and area of spherical segment between user bounds
  Real Vsh = 4./3. * PI * rsh*rsh*rsh * (cos(theta1)-cos(theta2));
  Real Ash = 4. * PI * rsh*rsh * (cos(theta1)-cos(theta2));
  // useful intermediate calculations
  Real rhoin = (E0 * 2*gam*gm1) / (t_stop*Ash * (2/Mach/Mach + gam*gm1)) / (vr*vr*vr);
  Real e_th = (E0 * 2) / (t_stop*Ash) / (2. + gam*gm1 * Mach*Mach) / vr; // calculate thermal energy density from mach number and total energy
  Real epsilon_k=(1./2. * gam*gm1 * Mach*Mach);

  std::cout << pcoord->x1f(is) << std::endl;
  std::cout << rsh << std::endl;
  // std::cout << "Initial conditions:  rhoin = " << rhoin << "  e_th = " << e_th << "  eps_k = " << epsilon_k << std::endl;
  // std::cout << "Initial conditions:  E0 = " << E0 << "  Ash = " << Ash << "  t_stop = " << t_stop << "  Mach = " << Mach << std::endl;
  // std::cout << "Initial conditions:  rsh = " << rsh << "  th1 = " << theta1 << "  th2 = " << theta2 << std::endl;

  // setup uniform ambient medium
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      Real theta_calc = pcoord->x2v(j);
      for (int i=is; i<=ie; i++) {
        // constant density medium
        phydro->u(IDN,k,j,i) = rhoa;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pa/gm1; //gamma/gm1;
        // putting all the energy between theta1 and theta2
        if (theta1 < theta_calc && theta_calc < theta2){
          // vr = sqrt(gam*pa*prat/rhoa/rhorat)*Mach;
          // vr = sqrt(gam)*Mach;
          phydro->u(IM1,k,j,is) = rhoin*vr;
          phydro->u(IDN,k,j,is) = rhoin;
          phydro->u(IEN,k,j,is) = e_th*(1.+epsilon_k);
        }
      }
    }
  }
}


void JetInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real gm1 = gam - 1.0;
  Real rsh = pco->x1v(il); // rs;  // this gives xi = 1.23
  Real Vsh = 4./3. * PI * rsh*rsh*rsh * (cos(theta1)-cos(theta2));
  Real Ash = 4. * PI * rsh*rsh * (cos(theta1)-cos(theta2));
  Real rhoin = (E0 * 2*gam*gm1) / (t_stop*Ash * (2/Mach/Mach + gam*gm1)) / (vr*vr*vr);
  Real e_th = (E0 * 2) / (t_stop*Ash) / (2. + gam*gm1 * Mach*Mach) / vr;
  Real pin = e_th*gm1;
  // std::cout << "Boundary injection:  rhoin = " << rhoin << "  e_th = " << e_th << "  pin = " << pin << std::endl;
  // std::cout << "Boundary injection:  E0 = " << E0 << "  Ash = " << Ash << "  t_stop = " << t_stop << "  Mach = " << Mach << std::endl;
  // std::cout << "Boundary injection:  rsh = " << rsh << "  th1 = " << theta1 << "  th2 = " << theta2 << std::endl;

  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      Real theta_calc = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
      //for (int i=0; i<=ngh; ++i) {
        // putting all the energy between theta1 and theta2
        // and spreading inside of radius rs
        if (theta1 < theta_calc && theta_calc < theta2){ // && pco->x1v(i) < rs){
          // vr = sqrt(gam*pa*prat/rhoa/rhorat)*Mach;
          if (time < t_stop){
            // translating from conseved to primitive variables
            pmb->phydro->w(IDN,k,j,il-i) = rhoin;
            pmb->phydro->w(IVX,k,j,il-i) = vr;
            pmb->phydro->w(IVY,k,j,il-i) = 0.0;
            pmb->phydro->w(IVZ,k,j,il-i) = 0.0;
            pmb->phydro->w(IPR,k,j,il-i) = pin;
          } else {
            pmb->phydro->w(IDN,k,j,il-i) = rhoin/100;
            pmb->phydro->w(IVX,k,j,il-i) = vr/100;
            pmb->phydro->w(IVY,k,j,il-i) = 0.0;
            pmb->phydro->w(IVZ,k,j,il-i) = 0.0;
            pmb->phydro->w(IPR,k,j,il-i) = pin/100;
          }

        } else {
          pmb->phydro->w(IDN,k,j,il-i) = pmb->phydro->w(IDN,k,j,il);
          pmb->phydro->w(IVX,k,j,il-i) = pmb->phydro->w(IVX,k,j,il);
          pmb->phydro->w(IVY,k,j,il-i) = pmb->phydro->w(IVY,k,j,il);
          pmb->phydro->w(IVZ,k,j,il-i) = pmb->phydro->w(IVZ,k,j,il);
          pmb->phydro->w(IPR,k,j,il-i) = pmb->phydro->w(IPR,k,j,il);
        }
      }
    }
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real time = pmy_mesh->time;
  std::cout << "[MeshBlock::UserWorkBeforeOutput]" << std::endl;
  gam = peos->GetGamma();
  Real gm1 = gam - 1.0;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is-NGHOST; i<=ie; i++) {

        Real ke = 0.5*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        Real e_gas = phydro->u(IEN,k,j,i) - ke;
        Real p = e_gas*gm1;

        // do something std::cout << rhoin << std::endl;
        // std::cout << "cons " << phydro->u(IDN,k,j,i) << ",  prim " << phydro->w(IDN,k,j,i) << std::endl;

        // do something std::cout << rhoin << std::endl
        // for (int idx=4; idx<=4; idx++)
        //   std::cout << "cons[" << idx << "] = " << phydro->u(idx,k,j,i) << ",  prim[" << idx << "] = " << phydro->w(idx,k,j,i) << ", p_byhand = " << p << std::endl;

        if (p < 0.)
            std::cout << "i = " << i << ", cons[" << IEN << "] = " << phydro->u(IEN,k,j,i) << ",  prim[" << IEN << "] = " << phydro->w(IEN,k,j,i) << ", p_byhand = " << p << std::endl;

      }
    }
  }
}

void zero_mom_i(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  for (int k = ks; k <=ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is-1; i >= is-ngh; --i) {

        prim(IDN,k,j,i) = prim(IDN,k,j,is);
        prim(IPR,k,j,i) = prim(IPR,k,j,is);
        prim(IM1,k,j,i) = 0.;
        prim(IM2,k,j,i) = 0.;
        prim(IM3,k,j,i) = 0.;
      }
    }
  }
  return;
}

void cgbc_i(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  Real     rs = pco->x1v(is);
  Real   rsp1 = pco->x1v(is+1);
  Real inv_dr = 1./(rsp1 - rs);

  Real  d_s,  d_sp1, d_slope;
  Real vp_s, vp_sp1, vp_slope;
  Real  p_s,  p_sp1, p_slope;
  Real bt_s, bt_sp1, bt_slope;
  Real bp_s, bp_sp1, bp_slope;

  for (int k = ks; k <=ke; ++k) {
    for (int j = js; j <= je; ++j) {

      d_s      = prim(IDN,k,j,is);
      d_sp1    = prim(IDN,k,j,is+1);

      p_s      = prim(IPR,k,j,is);
      p_sp1    = prim(IPR,k,j,is+1);

      vp_s     = prim(IM1,k,j,is);
      vp_sp1   = prim(IM1,k,j,is+1);

      d_slope  = inv_dr*( d_sp1 - d_s);
      vp_slope = inv_dr*(vp_sp1 - vp_s);
      p_slope  = inv_dr*( p_sp1 - p_s);

      if (MAGNETIC_FIELDS_ENABLED) {
        bt_s     = b.x2f(k,j,is);
        bt_sp1   = b.x2f(k,j,is+1);
        bp_s     = b.x3f(k,j,is);
        bp_sp1   = b.x3f(k,j,is+1);
        bt_slope = inv_dr*(bt_sp1 - bt_s);
        bp_slope = inv_dr*(bp_sp1 - bp_s);
      }

      for (int i = is-1; i >= is-ngh; --i) {

        Real r = pco->x1v(i);
        prim(IDN,k,j,i) = d_s + d_slope*(rs - r);
        prim(IPR,k,j,i) = p_s + p_slope*(rs - r);
        prim(IM1,k,j,i) = prim(IM1,k,j,is);
        prim(IM2,k,j,i) = prim(IM2,k,j,is);
        prim(IM3,k,j,i) = prim(IM3,k,j,is);

        if (MAGNETIC_FIELDS_ENABLED) {
          prim(IM3,k,j,i) = vp_s + vp_slope*(rs - r);
          b.x1f(k,j,i-1)  = b.x1f(k,j,i); // let the code calc B_r
          b.x2f(k,j,i)    = bt_s + bt_slope*(rs - r);
          b.x3f(k,j,i)    = bp_s + bp_slope*(rs - r);
        }
      }
    }
  }
  return;
}

Real shockProperties(MeshBlock *pmb, int iout) {

  int is=pmb->is;
  int ie=pmb->ie;
  int js=pmb->js;
  int je=pmb->je;
  int ks=pmb->ks;
  int ke=pmb->ke;

  Real    pMax = 0;
  Real  rhoMax = 0;
  int  pMaxLoc_r = 0;
  int  pMaxLoc_theta = 0;

  // bool startSum = false;
  bool startLocFound = false;
  int  EStartLoc = ie;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=ie; i>=is; --i) {
        Real p = pmb->phydro->w(IEN,k,j,i);
        if(p > pMax) {
          pMax = p;
          pMaxLoc_r = i;
          pMaxLoc_theta = j;
          rhoMax  = pmb->phydro->u(IDN,k,j,i);
        }
        if(!startLocFound && p>pa){
          startLocFound = true;
          EStartLoc     = i;
        }
      }
    }
  }

  if (iout == 1) { // return shock location in r
    return pmb->pcoord->x1v(pMaxLoc_r);
  }
  if (iout == 2) { // return shock location in theta
    return pmb->pcoord->x2v(pMaxLoc_theta);
  }
  if (iout == 3) { // return pressure max
    return pMax;
  }
  if (iout == 4) { // return density max
    return rhoMax;
  }
  if (iout == 5) { // return density max
    return pmb->phydro->w(IM1,ks,js,pMaxLoc_r);
  }

  // by default, return the energy of shocked region
  Real  Esh = 0.;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=EStartLoc; i>=is; --i) {
        Real  en = pmb->phydro->u(IEN,k,j,i);
        // Real  dV = pmb->pcoord->GetCellVolume(k,j,i);
        Real rout = pmb->pcoord->x1f(i+1);
        Real  rin = pmb->pcoord->x1f(i);
        Real   dV = 4./3. * PI * (rout*rout*rout - rin*rin*rin);
             Esh += en*dV;
      }
    }
  }
  return Esh;
}
