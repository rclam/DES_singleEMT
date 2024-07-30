#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic-read-temp.hpp"
#include "ic.hpp"

#include <cmath>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };

} // anonymous namespace


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return; // TAKE OUT COMMENT AFTER PURE SHEAR RUN
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    // ======================================================
    // Single EMT correction (one time only)
    double bulkm = var.mat->bulkm(0);
    double shearm = var.mat->shearm(0);
    double emt_rho = var.mat->emt_rho(0);
    double emt_pf = var.mat->emt_pf_z(0);
    double theta_normal = var.mat->theta_normal(0);

    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double E0 = shearm * ((3*lambda + 2*shearm)/(lambda + shearm)); // Young's Mod
    double v = lambda / (2*(lambda + shearm));                      // poisson ratio

    // ============ intact compliance S_i ============
    MatrixXd S_i{
        {1/E0, -v/E0, -v/E0, 0.0, 0.0, 0.0},
        {-v/E0, 1/E0, -v/E0, 0.0, 0.0, 0.0},
        {-v/E0, -v/E0, 1/E0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1/shearm, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1/shearm, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1/shearm}
    };
    // ============ Crack parameters ============
    double th_rad = ((90-theta_normal)+90)*(M_PI/180); // degree in radian, Uses cmath PI
    double a_n[3] = {cos(th_rad), sin(th_rad), 0.0};
    // crack density tensor alpha
    double a_alpha[3][3] = {
        {0.0,0.0,0.0},
        {0.0,0.0,0.0},
        {0.0,0.0,0.0}
        };
    // ============ Solving Alpha Tensor ============
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            a_alpha[i][j] = emt_rho*a_n[i]*a_n[j]; //tensor product of crack normal tensor
        }
    }
    // ============ Compliance Correction term  ============
    // correction term (S_voigt): delta_s_a * delta_s_b
    double delta_s_a = (8.0*(1.0-pow(v,2)))/(3.0*E0*(2.0-v));
    MatrixXd delta_s_b{
        {4*a_alpha[0][0], 0.0, 0.0, 0.0, 2*a_alpha[0][2], 2*a_alpha[0][1]}, 
        {0.0, 4*a_alpha[1][1], 0.0, 2*a_alpha[1][2], 0.0, 2*a_alpha[1][0]}, 
        {0.0, 0.0, 4*a_alpha[2][2], 2*a_alpha[2][1], 2*a_alpha[2][0], 0.0}, 
        {0.0, 2*a_alpha[2][1], 2*a_alpha[1][2], a_alpha[1][1] + a_alpha[2][2], a_alpha[1][0], a_alpha[2][0]}, 
        {2*a_alpha[2][0], 0.0, 2*a_alpha[0][2], a_alpha[0][1], a_alpha[0][0] + a_alpha[2][2], a_alpha[2][1]}, 
        {2*a_alpha[1][0], 2*a_alpha[0][1], 0.0, a_alpha[0][2], a_alpha[1][2], a_alpha[0][0]+a_alpha[1][1]}
    };
    //correction term
    MatrixXd S_voigt(6,6);
    S_voigt << delta_s_a * delta_s_b;

    // ============ Solving S_e ============
    // New Cracked Compliance S_e
    MatrixXd S_e(6,6);
    S_e << S_i + S_voigt;

    // ============ Solving c_e ============
    // New Cracked Stiffness c_e
    MatrixXd c_e(6,6);
    c_e << S_e.inverse();

    MatrixXd S_klmm{
        {S_i(0,0)+S_i(0,1)+S_i(0,2), S_i(5,0)+S_i(5,1)+S_i(5,2), S_i(4,0)+S_i(4,1)+S_i(4,2)},
        {S_i(5,0)+S_i(5,1)+S_i(5,2), S_i(1,0)+S_i(1,1)+S_i(1,2), S_i(3,0)+S_i(3,1)+S_i(3,2)},
        {S_i(4,0)+S_i(4,1)+S_i(4,2), S_i(3,0)+S_i(3,1)+S_i(3,2), S_i(2,0)+S_i(2,1)+S_i(2,2)}
    };

    VectorXd S_klmm_V(6);
    S_klmm_V(0) = S_klmm(0,0);
    S_klmm_V(1) = S_klmm(1,1);
    S_klmm_V(2) = S_klmm(2,2);
    S_klmm_V(3) = S_klmm(1,2);
    S_klmm_V(4) = S_klmm(0,2);
    S_klmm_V(5) = S_klmm(0,1);

    // ============ complianceXstiffness (product) ============
    // dot product CS in voigt
    VectorXd CS_V(6);
    CS_V = c_e*S_klmm_V;

    // CS voigt to full
    MatrixXd CS{
        {CS_V(0), CS_V(5), CS_V(4)},
        {CS_V(5), CS_V(1), CS_V(3)},
        {CS_V(4), CS_V(3), CS_V(2)}
    };
    // ============ Biot calc ============
    // B_ij = Kronecker - CS
    MatrixXd Kronecker_delta{
        {1.0,0.0,0.0},
        {0.0,1.0,0.0},
        {0.0,0.0,1.0}
    };

    MatrixXd Biot(3,3);
    //Biot << Kronecker_delta - CS;
    if(emt_rho == 0)
        Biot << Kronecker_delta;
    else
        Biot << Kronecker_delta - CS;
    
    // ============ stress_corr[i][j] ============
    MatrixXd stress_corr(3,3);
    stress_corr = emt_pf * Biot;
    //stress_corr = emt_pf * Kronecker_delta; //use to make effective elastoplastic run

    //       [s0 s3 s4]
    //   S = [s3 s1 s5]
    //       [s4 s5 s2]
    
    double stress_corr_V[NSTR] = {0.0};
    for (int i = 0; i < NDIMS; i++)
        stress_corr_V[i] = stress_corr(i,i);
    if( NDIMS == 2)
        stress_corr_V[2] = stress_corr(0,1);
    else if( NDIMS == 3) {
        stress_corr_V[5] = stress_corr(1,2);
        stress_corr_V[4] = stress_corr(0,2);
        stress_corr_V[3] = stress_corr(0,1);
    }
    // END of EMT correction calculation
    // ======================================================
    
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        double p = ref_pressure(param, zcenter);
        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = -p;
            strain[e][i] = -p / ks / NDIMS;
            //std::cerr << "\nstress[e][i]: " << stress[e][i];
            //std::cerr << "\nstrain[e][i]: " << strain[e][i];
        }
        if (param.mat.is_plane_strain)
            stressyy[e] = -p;
    }

    /*for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        double p = ref_pressure(param, zcenter);
        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        
        stress[e][0] = -105e6;
        strain[e][0] = -105e6 / ks / NDIMS;
        
        stress[e][1] = -25e6;
        strain[e][1] = -25e6 / ks / NDIMS;
        
        if (param.mat.is_plane_strain)
            stressyy[e] = -25e6;
    }*/

    if (var.mat->rheol_type & MatProps::rh_emt){
        std::cerr << "\nUpdating initial stress state for EMT\n";
        std::cerr << "stress_corr_V[0]: " << stress_corr_V[0];
        std::cerr << "\nstress_corr_V[1]: " << stress_corr_V[1];
        std::cerr << "\nstress_corr_V[2]: " << stress_corr_V[2] << "\n\n";
        for (int e=0; e<var.nelem; ++e) {
            // EMT correction diagonal stress
            for (int i=0; i<NDIMS; ++i){
                stress[e][i] += stress_corr_V[i]; // += is convention from rheology.cxx
                //std::cerr << "\nstress[e][i]: " << stress[e][i];
            }
            // EMT correction off-diagonal stress        
            for (int i=NDIMS; i<NSTR; ++i){
                stress[e][i] += stress_corr_V[i];
                //std::cerr << "\nstress_xy: " << stress[e][i];
            }
        }
    }

    compensation_pressure = ref_pressure(param, -param.mesh.zlength);
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;

    // TODO: adding different types of weak zone
    double plane_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        break;
    case 1:
        // a planar weak zone, cut through top center
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(plane_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(plane_center, semi_axis);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain;

        // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    delete weakzone;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    switch(param.ic.temperature_option) {
    case 0:
        {
            const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

            for (int i=0; i<var.nnode; ++i) {
                double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * age);
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }
    case 90:
        read_external_temperature_from_comsol(param, var, *var.temperature);
        break;
    default:
        std::cout << "Error: unknown ic.temperature option: " << param.ic.temperature_option << '\n';
        std::exit(1);
    }
}


