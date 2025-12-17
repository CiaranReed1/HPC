#pragma once
#define _XOPEN_SOURCE

const int M = 1024;  // x domain size
const int N = 512;   // y domain size
const int T = 50000; // number of time steps
const double uhi = 0.5;
const double ulo = -0.5;
const double vhi = 0.1;
const double vlo = -0.1;
const double a = 0.3;  // model parameter a
const double b = 0.1;  // model parameter b
const double c = 0.01; // model parameter c
const double d = 0.0;  // model parameter d
const double Ic = 0.03;
const double R = 1.0;
const double dt = 0.0625;          // time step
const double dx = 2.0;             // spatial resolution
const double DD = 1.0 / (dx * dx); // diffusion scaling
const int m = 100;                 // Norm calculation period

void init(double *u, double *v);

void step(const double *du, const double *dv, double *u, double *v);

void dxdt(double *du, double *dv, const double *u, const double *v);

double norm(const double *x);

double stim(int i, int j) {
  if ((i < 10) && (j < 10)) {
    return Ic;
  } else {
    return 0.0;
  }
}

inline double f(double u, double v) { return u * (1.0 - u) * (u - b) - v; }

inline double g(double u, double v) { return c * (a * u - v); }
