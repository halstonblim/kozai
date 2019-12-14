#include "kozai.h"

const string help="\nCarl's Kozai Code\n\n\
For integrating the secular equations of motions\n\
for the three-body problem, including the quadrupole\n\
and octupole terms (see Tremaine+2009 and Liu+2015)\n\n\
Also includes relativistic terms for the inner binary, including:\n\
  Pericenter precession (ala Ford+2000)\n\
  Gravitational Wave emission (ala Blaes+2002)\n\
  Spin-Orbit coupling (ala me, following Barker and O'Connell 1975)\n\
  Spin-Spin coupling (same)\n\n\
Options:\n\
  --no-int  don't integrate the system, just print info\n\
  --help    print this\n\n\
Parameters:\n\
  --m1      mass 1 of inner binary (MSUN)\n\
  --m2      mass 2 of inner binary (MSUN)\n\
  --m3      mass of tertiary (MSUN)\n\
  --a1      Semi-major axis of inner binary (AU)\n\
  --a2      Semi-major axis of outer binary (AU)\n\
  --e1      Eccentricity of inner binary\n\
  --e2      Eccentricity of outer binary\n\
  --g1      Argument of pericenter for inner binary (deg)\n\
  --g2      Argument of pericenter for outer binary (deg)\n\
  --omega1  Longitude of ascending node for inner binary (deg)\n\
  --omega2  Longitude of ascending node for outer binary (deg)\n\
  --inc     Mutual Inclination (deg)\n\
  --rad1    Radius of mass 1 (RSUN)\n\
  --rad2    Radius of mass 2 (RSUN)\n\
  --time    Time to integrate (years)\n\
  --dt      How often to print output\n\
                Positive dt outputs every dt years\n\
                Negative dt outputs every 10^{-dt} adaptive timesteps\n\
                    (e.g. dt=-2 outputs every 100 timesteps)\n\n\
Flags:\n\
  --quad       Include quadrupole terms \n\
  --oct        Include octupole Terms\n\
  --peri       Include inner pericenter precession (1pN) terms\n\
  --outerperi  Include outer pericenter precession (1pN) terms\n\
  --ignore_gsl Ignore errors from GSL (i.e. when a binary merges, the inner\n\
               binary can decouple from the outer s.t. the optimal timestep\n\
               for the inner is beyond machine-tolerance away from the outer)\n\
  --print_gsl  Print status to cour and state to cerr when GSL_ERROR flagged\n";

// Function for printing the state of a triple; placed here for easy
// modification
void print_state(double t_next, kozai_struct *kozai, ostream& stream = cerr){

  streamsize ss = stream.precision();

  stream << setprecision(14) << (t_next/YEAR) << "  " 
    << setprecision(8) 
    << kozai->get_a1()/AU << "  "
    << kozai->get_a2()/AU << "  " 
    << kozai->get_ecc1() << "  "
    << kozai->get_ecc2() << "  " 
    << setprecision(14) 
    << kozai->get_inc()/DEG << "  " 
    << kozai->get_inc1()/DEG << "  " 
    << kozai->get_inc2()/DEG << "  " 
    << setprecision(8) 
    << kozai->get_g1()/DEG << "  " 
    << kozai->get_g2()/DEG << "  "
    << kozai->get_Omega1()/DEG << "  "
    << kozai->get_Omega2()/DEG << endl;
}

int main(int argc, char **argv){

  double eps_abs = 1.e-12;
  double eps_rel = 1.e-12;
  double h = 1.e-10;
  double tmin=0, tmax=3.2e8*YEAR, delta_t=1e5*YEAR;
  bool collision=false;
  bool failure=false;
  bool IGNORE_GSL_ERRORS=false;
  bool PRINT_GSL=false;
  bool TENhz_ecc=false;
  bool TENhz_circ=false;
  double t=tmin;
  double t_insp;
  double max_e=0.;
  int status;

  //First, initialize the Kozai structure
  kozai_struct *kozai = new kozai_struct;
  set_parameters(argc, argv, kozai, tmax, delta_t, IGNORE_GSL_ERRORS, PRINT_GSL);
  kozai->initialize();

  //Then, set up the GSL integrator (we're using an 8th-order Runge Kutta)
  const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc(type_ptr, DIMENSION);
  gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc(DIMENSION);
  gsl_odeiv2_system my_system;

  my_system.function = rhs;
  my_system.dimension = DIMENSION;
  my_system.params = kozai;

  print_header_and_initial_state(kozai);
  print_state(0.,kozai);


  //Then actually integrate the system forward over each timestep
  long step = 0;
  long delta_output = 1;
  bool output_every_nsteps = false;

  for (double t_next = delta_t; t_next < tmax; t_next += delta_t){

    // if printing every x timesteps, the while loop does the whole
    // integration
    if (delta_t <= 0){
      delta_output = pow(10,-delta_t);
      delta_t = 1e100;
      t_next = tmax;
      output_every_nsteps = true;
    }

    //Then integrate the system forward for this timestep
    while (t < t_next){
      status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                   &my_system, &t, t_next, &h, kozai->get_y());

      //Check for collisions
      if(kozai->collision()){
        collision = true;
        break;
      }

      //Record the maximum ecc of the inner binary
      if (kozai->get_ecc1() > max_e)
        max_e = kozai->get_ecc1();
      
      //check for hitting the LIGO band
      if(kozai->gwave_freq() > 10 && TENhz_ecc == false){
        TENhz_ecc = true;
        cout << "LIGO band (eccentric): ";
        print_state(t,kozai,cout);
      }
      
      if(2*kozai->get_forb() > 10){
        TENhz_circ = true;
        cout << "LIGO band (circular): ";
        print_state(t,kozai,cout);
        break;
      }

      //print status if PRINT_GSL flagged
      if ((status != GSL_SUCCESS && PRINT_GSL) ){
        cout << "GSL_ERROR: (time, status) = " << t / YEAR << "," << gsl_strerror(status) << endl;
        print_state(t,kozai);
      }
      //check for some type of failure/NaN
      if ((status != GSL_SUCCESS && !IGNORE_GSL_ERRORS) || isnan(kozai->get_ecc1()) || kozai->get_ecc1() >= 1.){
        failure = true;
        break;
      }

      //otherwise, print the step (if we're printing on the integrator
      //timescale)
      if (output_every_nsteps && step % delta_output == 0)
        print_state(t,kozai);
        /*
        cout << "t = " << t/YEAR << endl;
        for(int z = 0;z <= 12; z++) cout << kozai->get_y()[z] << " ";
        cout << endl;
        for(int z = 0;z <= 12; z++) cout << evolve_ptr->yerr[z] << " " ;
        cout << endl << endl;
        */
      step += 1;
    }

    //Finally, print the state of the system
    print_state(t,kozai);

    if (collision == true){
      cout << "Collision!" << endl;
      break;
    }

    if (TENhz_circ == true){
      cout << "LIGO band!" << endl;
      break;
    }

    if (failure == true){
      cout << "NaNs (if near merger, try --ignore_gsl)!" << endl;
      break;
    }
  }

  streamsize ss = cout.precision();
  cout << setprecision(10);
  cout << "Maximum Eccentricity: " << max_e << endl;
  cout << setprecision(ss);

  /* all done; free up the gsl_odeiv stuff */
  gsl_odeiv2_evolve_free (evolve_ptr);
  gsl_odeiv2_control_free (control_ptr);
  gsl_odeiv2_step_free (step_ptr);
  delete kozai;
}

void print_header_and_initial_state(kozai_struct* kozai){
  cout << "Triple setup:\n"  
    << "  m1 = " << kozai->get_m1()/MSUN << endl
    << "  m2 = " << kozai->get_m2()/MSUN << endl
    << "  m3 = " << kozai->get_m3()/MSUN << endl
    << "  a1 = " << kozai->get_a1()/AU << endl
    << "  a2 = " << kozai->get_a2()/AU << endl
    << "  e1 = " << kozai->get_ecc1() << endl
    << "  e2 = " << kozai->get_ecc2() << endl
    << "  g1 = " << kozai->get_g1()/DEG << endl
    << "  g2 = " << kozai->get_g2()/DEG << endl
    << "  Omega1 = " << kozai->get_Omega1()/DEG << endl
    << "  Omega2 = " << kozai->get_Omega2()/DEG << endl
    << "  inc = " << kozai->get_inc()/DEG << endl
    << "  inc 1 = " << kozai->get_inc1()/DEG << endl
    << "  inc 2 = " << kozai->get_inc2()/DEG << endl
    << "  radius 1 = " << kozai->get_r1()/RSUN << endl
    << "  radius 2 = " << kozai->get_r2()/RSUN << endl;
  
  //Print the header with units and column names, and the initial state of the
  //system
  // cerr << "#1:t[yr]  #2:a1[AU]  #3:a2[AU] #4:e1  #5:e2  #6:inc[deg] #7:i1  #8:i2  #9:g1  #10:g2" << endl;
  cerr << "#1:t[yr]  #2:a1[AU]  #3:a2[AU] #4:e1  #5:e2  #6:inc[deg] #7:i1  #8:i2  #9:g1  #10:g2  #11:h1  #12:h2" << endl;
}


int rhs(double t, const double y[], double f[], void *kozai_ptr){

  kozai_struct *kozai = (kozai_struct *) kozai_ptr;
  //First, convert the y array into vectors
  //
  //Note: you can't use the vectors in the kozai structure; for adaptive
  //timstepping, GSL sometimes uses a temporary y array, which won't be
  //updated in the kozai class
  vec j1 = vec(y[0],y[1],y[2]);
  vec e1 = vec(y[3],y[4],y[5]);
  vec j2 = vec(y[6],y[7],y[8]);
  vec e2 = vec(y[9],y[10],y[11]);

  //extract semi-major axes, masses and angular momenta
  double a1 = y[12];
  double a2 = kozai->get_a2();
  double m1 = kozai->get_m1();
  double m2 = kozai->get_m2();
  double m3 = kozai->get_m3();
  double L1 = kozai->get_L1_no_a()*sqrt(a1);
  double L2 = kozai->get_L2_no_a()*sqrt(a2);

  //pre-compute a bunch of the vectorial quantities that are used repeatedly
  double j1n = abs(j1);
  double j2n = abs(j2);
  vec n1 = j1/j1n;
  vec n2 = j2/j2n;

  double e1n = abs(e1);
  double e2n = abs(e2);
  vec u1 = e1/e1n;
  vec u2 = e2/e2n;

  double j1n2 = j1*n2;
  double j1u2 = j1*u2;
  double e1n2 = e1*n2;
  double e1u2 = e1*u2;

  vec j1xn2 = j1^n2;
  vec j1xu2 = j1^u2;
  vec e1xn2 = e1^n2;
  vec e1xu2 = e1^u2;

  double L1L2 = L1 / L2;

  //Compute the (local) quadrupole timescale
  double tsec = (sqrt((m1+m2)/(G*pow(a1,3.)))*pow(a2,3.)*pow(1.-sqr(e2n),1.5) /m3);

  vec dj1dt=0., de1dt=0., dj2dt=0., de2dt=0.;
  double dadt  = 0;

  //Add the quadrupole-order secular evolution equations
  if(kozai->get_quadrupole()){
    double quad_coef = 0.75 / tsec;
    dj1dt += quad_coef*(j1n2*j1xn2 - 5.*e1n2*e1xn2);
    de1dt += quad_coef*((j1n2*e1xn2) + (2.*(j1^e1)) - (5.*(e1n2*j1xn2)));
    dj2dt += quad_coef*L1L2*((-j1n2*j1xn2) + (5*e1n2*e1xn2));
    de2dt += quad_coef*(L1L2/j2n)*((j1n2*(e2^j1)) - (5.*e1n2*(e2^e1)) 
         - ((0.5 - 3.*sqr(e1n) + 12.5*sqr(e1n2) - 2.5*sqr(j1n2))*(n2^e2)));
  }

  //Add the octupole-order secular evolution equations
  if(kozai->get_octupole() == true){
    double octo_coef = -1.171875*(fabs(m1-m2)/(m1+m2)*(a1/a2)*(e2n/sqr(j2n)))/tsec;
    dj1dt += octo_coef*((((2.*(e1u2*j1n2+e1n2*j1u2)*j1)
        + 2.*(j1u2*j1n2-7.*e1u2*e1n2)*e1)^n2) + ((((2.*e1n2*j1n2)*j1)
        + (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e1)^u2));
    de1dt += octo_coef*(((((2.*e1n2*j1n2)*e1) + (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)   
        + sqr(j1n2))*j1)^u2) + (((2.*(e1u2*j1n2+e1n2*j1u2)*e1) 
        + (2.*(j1n2*j1u2-7*e1n2*e1u2)*j1))^n2) + (3.2*e1u2*(j1^e1)));
    dj2dt += octo_coef*L1L2*((2*(((e1n2*j1u2)*n2) + ((e1u2*j1n2)*n2)
        + (e1n2*j1n2)*u2)^j1) + ((((2*j1u2*j1n2)*n2) - ((14.*e1u2*e1n2)*n2) 
        + ((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*u2))^e1));
    de2dt += octo_coef*(L1L2/j2n)*(2*(((e1n2*((j1*e2)*u2))
        + (j1n2*((e1*e2)*u2)) + (((sqr(j2n)/e2n)*e1n2*j1n2)*n2))^j1) + ((((2.*e2n*j1u2*j1n2)*u2)
        - ((14.*e1n2*(e1*e2))*u2) + (sqr(j2n)/e2n)*((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)
        + sqr(j1n2))*n2))^e1) - (((2.*(0.2-1.6*sqr(e1n))*e1u2*e2) + ((14.*e1n2*j1u2*j1n2)*e2)
        + (7.*e1u2*(1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e2))^n2));
  }

  //Add the pericenter precession of the inner binary
  if (kozai->get_pericenter() == true)
    de1dt += (3./(c*c*a1*sqr(j1n)))*(pow(G*(m1+m2)/a1, 1.5) * (n1^e1));

  //Add the pericenter precession of the outer binary
  if (kozai->get_outer_pericenter() == true)
    de2dt += (3./(c*c*a2*sqr(j2n)))*(pow(G*(m1+m2+m3)/a2, 1.5) * (n2^e2));

  //Add 1pn-quadrupole cross terms
  if (kozai->get_cross_lim() == true){   
    // precompute various quantities
    vec v1 = n1^u1;
    vec v2 = n2^u2;
    double p2    = a2 * sqr(j2n);
    double m = m1+m2;
    double eta = m1*m2/sqr(m);
    double mtot = m + m3;
        
    // Define trig funcs wrt to inc, i1
    // Avoid evaluating angles; use vectors; Naoz et al (2013a) Appendix A
    double cinc    = n1*n2;
    double sinc    = sqrt(1.-sqr(cinc));
    double s2inc   = 2. * sinc * cinc;
    double sincsq  = sqr(sinc);
    double cincsq  = sqr(cinc);
    double c2inc   = 1 - 2*sqr(sinc);
    double cinc1   = n1[2];
    // double cinc2   = n2[2];
    double sinc1   = sqrt(1.-sqr(cinc1));
    double cscinc1 = 1 / sinc1;
    // double sinc2   = sqrt(1.-sqr(cinc2));
    double cotinc1     = 1 / sinc1;
    double cotinc1sinc = cinc1 * sinc / sinc1;
    double cscinc1sinc = sinc / sinc1;

    // Define trig functions wrt to g1, g2; Goes bad as sinc -> 0
    double cg1 = (n2*v1) / sinc;
    double sg1 = (n2*u1) / sinc;
    double cg2 = (n1*v2) / sinc;
    double sg2 = (n1*u2) / sinc;
    double c2g1 = 1 - 2*sqr(sg1);
    double s2g1 = 2 * sg1 * cg1;
    double c2g2 = 1 - 2*sqr(sg2);
    double s2g2 = 2 * sg2 * cg2;

    double cinc2 = n2[2];
    double sinc2 = sqrt(1.-sqr(cinc2));
    // 1PN-quadrupole cross terms; Lim Rodriguez (2019)
    // Inner binary
    double de1dtcross=0, dp1dtcross=0, dg1dtcross=0, dh1dtcross=0, di1dtcross=0;    
     
    // (nr, nm) = (1.5, -0.5) 
    
    double factor = pow(a1,1.5)*pow(a2,-4)*pow(c,-2)*pow(G,1.5)*pow(j2n,-5)*pow(m,-0.5)*pow(mtot,2);
    de1dtcross += factor * (-15*e1n*j1n*(8*c2g1*cinc*s2g2 - 2*s2g1*(c2g2*(3 + c2inc) + 6*sincsq))*pow(e2n,2))/32.;
    dp1dtcross += factor * (-15*a1*j1n*(-4*c2g1*cinc*s2g2 + s2g1*(c2g2*(3 + c2inc) + 6*sincsq))*pow(e1n,2)*pow(e2n,2))/8.;
    di1dtcross += factor * (-3*sinc*pow(e2n,2)*pow(j1n,-1)*(-5*(-3 + c2g2)*cinc*s2g1*pow(e1n,2) + s2g2*(5 + 5*c2g1*pow(e1n,2) - 3*pow(j1n,2))))/8.;
    dg1dtcross += factor * (-3*pow(j1n,-1)*(10*cotinc1sinc*s2g1*s2g2*pow(e1n,2)*(-16 + 3*pow(j2n,2)) + 20*cinc*s2g1*s2g2*pow(j1n,2)*(-16 + 3*pow(j2n,2)) - cotinc1*s2inc*(5 - 5*c2g1*pow(e1n,2) - 3*pow(j1n,2))*(4 - 10*pow(j2n,2) + c2g2*(-16 + 3*pow(j2n,2))) + (-3 + 5*c2g1)*c2inc*pow(j1n,2)*(4 - 10*pow(j2n,2) + c2g2*(-16 + 3*pow(j2n,2))) + (1 + 5*c2g1)*pow(j1n,2)*(-4 + 10*pow(j2n,2) + c2g2*(-48 + 9*pow(j2n,2)))))/64.;
    dh1dtcross += factor * (3*cscinc1*sinc*pow(e2n,2)*pow(j1n,-1)*(-5*s2g1*s2g2*pow(e1n,2) + (-3 + c2g2)*cinc*(5 - 5*c2g1*pow(e1n,2) - 3*pow(j1n,2)))*pow(j2n,-5))/8.;
    
    // (nr, nm) = (0, 0)
    dh1dtcross += (3*cscinc1sinc*pow(a2,-2.5)*pow(c,-2)*pow(G,1.5)*pow(j2n,-2)*pow(mtot,1.5))/2.;
    dg1dtcross += (3*(cinc - cotinc1sinc)*pow(a2,-2.5)*pow(c,-2)*pow(G,1.5)*pow(j2n,-2)*pow(mtot,1.5))/2.;
     
    // (nr, nm) = (-1, 1)
    dg1dtcross+= (-15*m*(4*c2g1*c2g2*cinc + (3 + c2inc)*s2g1*s2g2)*pow(a1,-1)*pow(a2,-1.5)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j1n,-3)*pow(j2n,-3)*(1 + j2n - 2*pow(j2n,2))*pow(1 + j2n,-1)*pow(mtot,0.5))/16.;
    

    // (nr, nm) = (1/2, 1/2)
    // de1dtcross += (3*j1n*(-11 + j1n*(-22 + (-23 + 12*eta)*j1n))*m3*s2g1*sincsq*pow(a1,0.5)*pow(a2,-3)*pow(c,-2)*pow(e1n,-3)*pow(G,1.5)*pow(-1 + j1n,2)*pow(j2n,-3)*pow(m,0.5))/8.;
    // dp1dtcross += (33*j1n*m3*s2g1*sincsq*pow(a1,1.5)*pow(a2,-3)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j2n,-3)*pow(m,0.5))/4.;
    // di1dtcross += (33*cg1*cinc*m3*sg1*sinc*pow(a1,0.5)*pow(a2,-3)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j1n,-1)*pow(j2n,-3)*pow(m,0.5))/4.;
    // dg1dtcross += (m3*pow(a1,0.5)*pow(a2,-3)*pow(c,-2)*pow(e1n,-4)*pow(G,1.5)*pow(-1 + j1n,2)*pow(j1n,-1)*pow(j2n,-3)*pow(m,0.5)*(3*c2g1*(-3 + c2inc)*(17 - 6*eta + j1n*(34 + 6*eta*(-2 + j1n) + 5*j1n))*pow(j1n,2) + 6*cotinc1*s2inc*(-11 + 11*c2g1*pow(e1n,2) + (5 - 4*eta)*pow(j1n,2))*pow(1 + j1n,2) + (-1 + 3*c2inc)*(-11 + 10*eta)*pow(j1n + pow(j1n,2),2) + 6*cincsq*(c2g1*(17 - 6*eta + j1n*(34 + 6*eta*(-2 + j1n) + 5*j1n))*pow(j1n,2) + (-11 + 10*eta)*pow(j1n + pow(j1n,2),2))))/32.;
    // dh1dtcross += (3*cscinc1*m3*s2inc*pow(a1,0.5)*pow(a2,-3)*pow(c,-2)*pow(G,1.5)*pow(j1n,-1)*(11 - 11*c2g1*pow(e1n,2) + (-5 + 4*eta)*pow(j1n,2))*pow(j2n,-3)*pow(m,0.5))/16.;

    de1dt += de1dtcross*u1 + (e1n*(dg1dtcross + dh1dtcross*cinc1))*v1 + (e1n*(di1dtcross*sg1-dh1dtcross*cg1*sinc1))*n1;
    dj1dt += (j1n*(-di1dtcross*sg1 + dh1dtcross*cg1*sinc1))*u1 + (-j1n*(di1dtcross*cg1 + dh1dtcross*sg1*sinc1))*v1 + ((-e1n/j1n)*de1dtcross)*n1;
    //de2dt += (e2n*(dh1dtcross*cinc2))*v2 + (e2n*(-1*dh1dtcross*cg2*sinc2))*n2;
    //dj2dt += (j2n*(dh1dtcross*cg2*sinc2))*u2 + (-j2n*(dh1dtcross*sg2*sinc2))*v2;
    dadt  += (dp1dtcross + 2. * a1 * e1n * de1dtcross) / sqr(j1n);

    // Outer binary terms are numerically insignificant
    // double de2dtcross=0, dp2dtcross=0, dg2dtcross=0, dh2dtcross=0, di2dtcross=0;    
    // // (nr, nm) = (2.0, -4.5)
    // de2dtcross += (3*eta*pow(a1,2)*pow(a2,-4.5)*pow(c,-2)*pow(e2n,-1)*pow(G,1.5)*(20*c2g2*cinc*s2g1*pow(e1n,2) - s2g2*(5*c2g1*(3 + c2inc)*pow(e1n,2) + 2*sincsq*(5 - 3*pow(j1n,2))))*pow(j2n,-6)*(25 - 7*pow(e2n,2)*(-2 + pow(j2n,2)) - 37*pow(j2n,2) + 12*pow(j2n,4))*pow(mtot,1.5))/128.;
    // di2dtcross += (3*eta*sinc*pow(a1,2)*pow(a2,-4.5)*pow(c,-2)*pow(G,1.5)*(5*(5 + 2*c2g2)*s2g1*pow(e1n,2) - 2*cinc*s2g2*(-5 + 5*c2g1*pow(e1n,2) + 3*pow(j1n,2)))*pow(j2n,-6)*(-1 + pow(j2n,2))*pow(mtot,1.5))/8.;
    // dp2dtcross += (-3*eta*pow(a1,2)*pow(a2,-3.5)*pow(c,-2)*pow(G,1.5)*(20*c2g2*cinc*s2g1*pow(e1n,2) - s2g2*(5*c2g1*(3 + c2inc)*pow(e1n,2) + 2*sincsq*(5 - 3*pow(j1n,2))))*pow(j2n,-4)*(-1 + pow(j2n,2))*pow(mtot,1.5))/4.;
    // dg2dtcross += (3*eta*pow(a1,2)*pow(a2,-4.5)*pow(c,-2)*pow(e2n,-2)*pow(G,1.5)*pow(j2n,-6)*pow(mtot,1.5)*(-16*cinc*cotinc2sinc*pow(cg1,2)*pow(cg2,2)*pow(e2n,2)*pow(j1n,2)*(-20 + pow(j2n,2)) - (-30*c2g1*sincsq*pow(e1n,2) + (1 + 3*c2inc)*(-5 + 3*pow(j1n,2)))*(-5 + pow(j2n,2))*(-4 + pow(j2n,2)) - 40*cotinc2sinc*s2g1*s2g2*pow(e1n,2)*pow(e2n,2)*(12 + pow(j2n,2)) - 40*cinc*s2g1*s2g2*pow(e1n,2)*(-44 + 40*pow(j2n,2) + pow(j2n,4)) - c2g2*(10*c2g1*(3 + c2inc)*pow(e1n,2) + 4*sincsq*(5 - 3*pow(j1n,2)))*(-44 + 40*pow(j2n,2) + pow(j2n,4)) + 80*cinc*s2g1*s2g2*pow(e1n,2)*(-32 + 33*pow(j2n,2) + 2*pow(j2n,4)) + 16*(-5 + 3*pow(j1n,2))*(-10 - 5*pow(j2n,2) + 3*pow(j2n,4)) + 8*cincsq*pow(cg1,2)*pow(cg2,2)*pow(j1n,2)*(34 - 81*pow(j2n,2) + 5*pow(j2n,4)) - 8*pow(cg1,2)*pow(cg2,2)*(-5 + 4*pow(j1n,2))*(-94 + 51*pow(j2n,2) + 13*pow(j2n,4)) + 16*cinc*cotinc2sinc*pow(cg2,2)*pow(e2n,2)*(-5 + 4*pow(j1n,2))*(-20 + pow(j2n,2))*pow(sg1,2) - 8*cincsq*pow(cg2,2)*(-5 + 4*pow(j1n,2))*(34 - 81*pow(j2n,2) + 5*pow(j2n,4))*pow(sg1,2) + 8*pow(cg2,2)*pow(j1n,2)*(-94 + 51*pow(j2n,2) + 13*pow(j2n,4))*pow(sg1,2) - 16*cinc*cotinc2sinc*pow(cg1,2)*pow(e2n,2)*pow(j1n,2)*(4 + 3*pow(j2n,2))*pow(sg2,2) - 8*pow(cg1,2)*(-5 + 4*pow(j1n,2))*(34 - 81*pow(j2n,2) + 5*pow(j2n,4))*pow(sg2,2) + 8*cincsq*pow(cg1,2)*pow(j1n,2)*(-94 + 51*pow(j2n,2) + 13*pow(j2n,4))*pow(sg2,2) + 16*cinc*cotinc2sinc*pow(e2n,2)*(-5 + 4*pow(j1n,2))*(4 + 3*pow(j2n,2))*pow(sg1,2)*pow(sg2,2) + 8*pow(j1n,2)*(34 - 81*pow(j2n,2) + 5*pow(j2n,4))*pow(sg1,2)*pow(sg2,2) - 8*cincsq*(-5 + 4*pow(j1n,2))*(-94 + 51*pow(j2n,2) + 13*pow(j2n,4))*pow(sg1,2)*pow(sg2,2) + 4*(2*(-5 + 3*pow(j1n,2))*(122 - 137*pow(j2n,2) + 27*pow(j2n,4)) + 40*cg1*cg2*sg1*sg2*pow(e1n,2)*(cotinc2sinc*pow(e2n,2)*(20 - 7*pow(j2n,2)) + 3*cinc*(6 - 14*pow(j2n,2) + 5*pow(j2n,4))) + pow(sg1,2)*(pow(cg2,2)*(cinc*(-5 + 4*pow(j1n,2))*(4*cotinc2sinc*pow(e2n,2)*(8 + 11*pow(j2n,2)) + cinc*(-330 + 327*pow(j2n,2) - 51*pow(j2n,4))) + 3*pow(j1n,2)*(134 - 165*pow(j2n,2) + 37*pow(j2n,4))) + (3*pow(j1n,2)*(110 - 109*pow(j2n,2) + 17*pow(j2n,4)) + cinc*(-5 + 4*pow(j1n,2))*(4*cotinc2sinc*pow(e2n,2)*(-32 + 25*pow(j2n,2)) - 3*cinc*(134 - 165*pow(j2n,2) + 37*pow(j2n,4))))*pow(sg2,2)) + pow(cg1,2)*(pow(cg2,2)*(-3*(-5 + 4*pow(j1n,2))*(134 - 165*pow(j2n,2) + 37*pow(j2n,4)) + cinc*pow(j1n,2)*(-4*cotinc2sinc*pow(e2n,2)*(8 + 11*pow(j2n,2)) + cinc*(330 - 327*pow(j2n,2) + 51*pow(j2n,4)))) + (-3*(-5 + 4*pow(j1n,2))*(110 - 109*pow(j2n,2) + 17*pow(j2n,4)) + cinc*pow(j1n,2)*(4*cotinc2sinc*pow(e2n,2)*(32 - 25*pow(j2n,2)) + 3*cinc*(134 - 165*pow(j2n,2) + 37*pow(j2n,4))))*pow(sg2,2)))))/256.;
    // dh2dtcross += (3*eta*cscinc2sinc*pow(a1,2)*pow(a2,-4.5)*pow(c,-2)*pow(G,1.5)*(10*s2g1*s2g2*pow(e1n,2) + (-5 + 2*c2g2)*cinc*(-5 + 5*c2g1*pow(e1n,2) + 3*pow(j1n,2)))*pow(j2n,-6)*(-1 + pow(j2n,2))*pow(mtot,1.5))/8.;

    // // (nr, nm) = (-0.5, 1.5)
    // dg2dtcross += (45*eta*(c2g1*c2g2*cinc + 2*cg1*cg2*(1 + cincsq)*sg1*sg2)*pow(a1,-0.5)*pow(a2,-2)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j1n,-2)*pow(j2n,-3)*pow(m,1.5))/16.;
    // de2dtcross += (-45*e2n*eta*(-4*c2g1*cinc*s2g2 + s2g1*(c2g2*(3 + c2inc) + 12*sincsq))*pow(a1,-0.5)*pow(a2,-2)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j1n,-2)*pow(j2n,-3)*pow(m,1.5))/64.;
    
    // de2dt += de2dtcross*u2 + (e2n*(dg2dtcross))*v2;
    // dj2dt += ((-e2n/j2n)*de2dtcross)*n2;
    // da2dt += (dp2dtcross + 2. * a2 * e2n * de2dtcross) / sqr(j2n);

    }

  //Add gravitational-wave emission for the inner binary
  if (kozai->get_radiation() == true){
    double dedt  = -(c304o15*e1n*(G3 * m1*m2*(m1+m2)) /
         (c5 * pow(a1,4.) * pow(j1n,5.))) * (1. + c121o304*sqr(e1n));
    dadt += -c64o5 * ((G3 * m1*m2*(m1+m2)) / (c5
        * pow(a1,3.) * pow(j1n,7.))) * (1. + c73o24*sqr(e1n) + c37o96*pow(e1n,4.)); 

    de1dt += dedt * u1;
    dj1dt += (-e1n / j1n) * dedt * n1;
  }

  //Finally, copy the derivatives into the output array
  for(int i=0; i<3; i++){
    f[i] = dj1dt[i];
    f[i+3] = de1dt[i];
    f[i+6] = dj2dt[i];
    f[i+9] = de2dt[i];
  }
  f[12] = dadt;

  return GSL_SUCCESS;
}

// Integrate the peters equation to find the inspiral time for a binary
// Not actually used at the moment, but probably useful to have around...
double peters_t(kozai_struct *kozai){
  double beta = (64./5.) * G3c5 * kozai->get_m1()*kozai->get_m2()*(kozai->get_m1()+kozai->get_m2());
  double e = kozai->get_ecc1();
  double c0 = kozai->get_a1() * (1-sqr(e)) * pow(e,-c12o19) * pow(1+(c121o304)*sqr(e),-c870o2299);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &peters_integral;
  double integral, error;

  gsl_integration_qags(&F,0,e,1e-7,1e-8,1000, w,&integral, &error);
  gsl_integration_workspace_free(w);

  return integral * c12o19 * sqr(c0)*sqr(c0) / beta;
}

double peters_integral(double e, void *params){
  return pow(e,c29o19) * pow(1+(c121o304)*sqr(e),c1181o2299) / pow(1-sqr(e),1.5);
}

void set_parameters(int argc, char **argv, kozai_struct *kozai, double &t_end, double &t_step, bool &IGNORE_GSL_ERRORS, bool &PRINT_GSL){
  int c;
  int opterr = 0;
  int option_index = 0;
  double radius;

  const struct option longopts[] =
  {
    {"m1",    1,  0, 'm'},
    {"m2",    1,  0, 'M'},
    {"m3",    1,  0, 'n'},
    {"a1",    1,  0, 'a'},
    {"a2",    1,  0, 'A'},
    {"e1",    1,  0, 'e'},
    {"e2",    1,  0, 'E'},
    {"g1",    1,  0, 'g'},
    {"g2",    1,  0, 'G'},
    {"omega1",1,  0, 'l'},
    {"omega2",1,  0, 'L'},
    {"inc"   ,1,  0, 'i'},
    {"rad1"   ,1,  0, 'b'},
    {"rad2"   ,1,  0, 'B'},
    {"quad"  ,0,  0, 'q'},
    {"oct"   ,0,  0, 'o'},
    {"peri"  ,0,  0, 'p'},
    {"outerperi",0,0,'P'},
    {"cross_lim",0,0,'X'},
    {"rad"   ,0,  0, 'r'},
    {"print_gsl",0,0,'f'},
    {"time"   ,1,  0, 'd'},
    {"dt"   ,1,  0, 'D'},
    {"ignore_gsl"   ,0,  0, 'I'},
    {"help"   ,0,  0, 'h'},
    {0,0,0,0},
  };

  while ((c = getopt_long (argc, argv, 
          "m:M:n:a:A:e:E:g:G:l:L:i:c:C:t:T:u:U:qopPXsSrfId:D:b:B:h",
          longopts,&option_index)) != -1)
  switch (c)
    {
    case 'm':
      kozai->set_m1(atof(optarg)*MSUN); break;
    case 'M':
      kozai->set_m2(atof(optarg)*MSUN); break;
    case 'n':
      kozai->set_m3(atof(optarg)*MSUN); break;
    case 'a':
      kozai->set_a1(atof(optarg)*AU); break;
    case 'A':
      kozai->set_a2(atof(optarg)*AU); break;
    case 'e':
      kozai->set_ecc1(atof(optarg)); break;
    case 'E':
      kozai->set_ecc2(atof(optarg)); break;
    case 'g':
      kozai->set_g1(atof(optarg)*DEG); break;
    case 'G':
      kozai->set_g2(atof(optarg)*DEG); break;
    case 'l':
      kozai->set_Omega1(atof(optarg)*DEG); break;
    case 'L':
      kozai->set_Omega2(atof(optarg)*DEG); break;
    case 'i':
      kozai->set_inc(atof(optarg)*DEG); break;
    case 'b':
      kozai->set_r1(atof(optarg)*RSUN); break;
    case 'B':
      kozai->set_r2(atof(optarg)*RSUN); break;
    case 'q':
      kozai->set_quadrupole(true); break;
    case 'o':
      kozai->set_octupole(true); break;
    case 'p':
      kozai->set_pericenter(true); break;
    case 'P':
      kozai->set_outer_pericenter(true); break;
    case 'X':
      kozai->set_cross_lim(true); break;
    case 'r':
      kozai->set_radiation(true); break;
    case 'f':
      PRINT_GSL = true; break;
    case 'I':
      IGNORE_GSL_ERRORS = true; break;
    case 'd':
      t_end = atof(optarg)*YEAR; break;
    case 'D':
      t_step = atof(optarg);
      if (t_step <= 0) break;
      t_step*=YEAR;
      break;
    case 'h':
      cout << help;
      exit(1);
    case '?':
      cout << help;
      exit(1);
    default:
      exit(1);
      }

}
