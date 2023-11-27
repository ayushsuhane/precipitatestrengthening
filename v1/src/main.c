/* Calculates the precipitation hardening behavior in Al-Mg system */
/*  Input 
        - heat treatment time, and temperature
        - alloy composition
        - solute diffusivity
        - solubility
        - shearable/non-shearable critical radius
    Outputs 
        - Strength for a given heat treatment
        - Maximum strength achieved during the heat treatment
*/

#include<stdio.h>
#include<math.h>

#define PI 3.14159

double calc_equilibrium_comp(double comp0, double deltaG_sol, 
                        double temperature, double gas_constant)
{
    return comp0 * exp(-deltaG_sol/(gas_constant * temperature));
}

double calc_diffusivity(double D0, double deltaG_d, 
                        double temperature, double gas_constant)
{
    return D0 * exp(-deltaG_d/(gas_constant * temperature));
}

double calc_drdt_nucleation(double radius, double N_nucleus, 
                                  double dndt, double alpha, 
                                  double r_star)
{
    return (1/N_nucleus)*dndt*(alpha*r_star - radius);
}

double calc_drdt_growth(double radius, double diffusivity, 
                        double comp, double comp_eq,
                        double r0)
{
    return (diffusivity/radius) * (comp - comp_eq * exp(r0/radius))/(1 - comp_eq * exp(r0/radius));
}


double calc_drdt_coarsening(double kappa, double radius)
{
    return kappa/(radius * radius);
}

double calc_dndt_growth(double zeldovich_parameter, double beta, 
                        double atomic_volume, double deltaG_star, 
                        double gas_constant, double temperature,
                        double tau, double time)
{
    return zeldovich_parameter * beta * atomic_volume * exp(- deltaG_star/(gas_constant * temperature)) * exp(-tau/time);
}

double calc_dndt_coarsening(double radius, double kappa, 
                            double r0, double comp, 
                            double N_nucleus)
{
    return (kappa/pow(radius, 3)) * ( (r0 * comp)/(radius * (1 - comp)) * ( 3.0 / (4 * PI * pow(radius, 3)) - N_nucleus) - 3 * N_nucleus);

}


double calc_fcoarse(double radius, double r_star)
{
    return 1 - erf(4.0 * (radius/r_star - 1));
}

double calc_transition(double dndt_growth, double dndt_coarsening)
{
    int flag;
    if (-dndt_coarsening > dndt_growth) flag = 1;
    else flag = 0;
    return flag;
}

double calc_volumefraction(double radius, double N_nucleus)
{
    return 4.0/3 * (PI * pow(radius, 3) * N_nucleus);
}

double calc_composition(double radius, double N_nucleus, double comp_initial)
{
    double vf = calc_volumefraction(radius, N_nucleus);
    return (comp_initial - vf)/(1 - vf);
}

int main(int argc, char *argv[])
{
    /*Input parameters*/
    double total_time = 1e6, temperature_inC = 160, comp_initial = 0.06;

    /*Internal Parameters*/
    double comp0 = 1.13, deltaG_sol = 17000; // equilibrium composition parameters
    double deltaG_nucleation = 168900; //eV
    double surface_energy = 0.3; // J/m^2
    double N_initial = 10000;
    double D0 = 2.348e-5, deltaG_d = 126211; //diffusion parameters, m^2/s, J/mol

    /*Variables for stress calculation*/
    double r_critical = 3.4e-9; //m
    double bulk_modulus = 26000; //MPa
    double burger_vector = 2.86e-10;
    double taylor_factor = 1.5;
    double sigma_0 = 10; //MPa
    double constant_solidsolution = 840; //MPa
    double line_force = 0.5 * bulk_modulus * burger_vector * burger_vector;



    /*Constants*/
    double gas_constant = 8.314;
    double lattice_parameter = 5.18e-10; // Angstroms
    double zeldovich_parameter = 0.05; 
    double molar_volume = 7.1e-6; //m^3/mol
    double Na = 6.023e23; // atoms/mol
    double atomic_volume = Na/molar_volume; // atoms/m^3


    /*Variables*/
    double time = 0;
    //double dt = 0.1;
    double dt_initial = 0.01;

    
    double alpha = 1.05;
    
    double N_nucleus = N_initial;
    double comp = comp_initial;

    double comp_new, r_new, N_new;
    double temperature = temperature_inC + 273;

    double diffusivity = calc_diffusivity(D0, deltaG_d, temperature, gas_constant);
    double comp_eq = calc_equilibrium_comp(comp0, deltaG_sol, temperature, gas_constant);
    double deltaG_vol, r0, r_star, radius, beta;
    r0 = 2.0 * surface_energy * molar_volume / (gas_constant * temperature);

    double kappa = (4.0/27) * comp_eq/(1-comp_eq) * r0 * diffusivity;
    double drdt_nucleation, drdt_growth, drdt_coarsening;
    double dndt_growth, dndt_coarsening;
    double fcoarse, transition;
    double drdt, dndt = 0;
    double vf, tau, deltaG_star, dt;

    double force, interparticle_spacing, sigma_ppt, sigma_solidsolution, sigma_total;
    double sigma_max = 0.0, time_max = 0.0;

    
    deltaG_vol = - gas_constant * temperature * log(comp/comp_eq); // J/mol
    r_star = -2.0 * surface_energy * molar_volume / deltaG_vol ;
    radius = r_star;

    int count = 0;


    while((time < total_time) && (comp > comp_eq) )
    {
        dt = pow(10, log10(dt_initial) + (count)/100.0);
        deltaG_vol = - gas_constant * temperature * log(comp/comp_eq); // J/mol
        r_star = -2.0 * surface_energy * molar_volume / deltaG_vol ;
        beta = 4 * PI * r_star * r_star * diffusivity * comp / pow(lattice_parameter, 4);
        deltaG_star = deltaG_nucleation/(log(comp/comp_eq) * log(comp/comp_eq)); // J/mol
        tau = 1.0/(2*beta*zeldovich_parameter);

    
        drdt_nucleation = calc_drdt_nucleation(radius, N_nucleus, dndt, alpha, r_star);
        drdt_growth = calc_drdt_growth(radius, diffusivity, comp, comp_eq, r0);
        drdt_coarsening = calc_drdt_coarsening(kappa, radius);
        
        dndt_growth = calc_dndt_growth(zeldovich_parameter, beta, atomic_volume, deltaG_star, gas_constant, temperature, tau, time);
        dndt_coarsening = calc_dndt_coarsening(radius, kappa, r0, comp, N_nucleus);
    
        fcoarse = calc_fcoarse(radius, r_star);
        transition = calc_transition(dndt_coarsening, dndt_growth);
        
        if(transition > 0.1) drdt = (1 - fcoarse)*drdt_growth + fcoarse*drdt_coarsening;
        else drdt = drdt_growth + drdt_nucleation;
        if(transition > 0.1) dndt = fcoarse*dndt_coarsening;
        else dndt = dndt_growth;
    
        r_new = radius + drdt*dt;
        N_new = N_nucleus + dndt*dt;
    
        comp_new = calc_composition(r_new, N_new, comp_initial);
        vf = calc_volumefraction(r_new, N_new);
        radius = r_new;
        N_nucleus = N_new;
        comp = comp_new;
        

        if(radius < r_critical) force = 2 * line_force * radius/r_critical;
        else force = 2 * line_force;

        interparticle_spacing = (radius / 2.0)  * sqrt(PI / vf);

        sigma_solidsolution = constant_solidsolution * pow(comp, 2.0/3.0);
        sigma_ppt = taylor_factor * pow(force/2.0, 1.5) / (burger_vector * pow(line_force, 0.5) * interparticle_spacing );
        sigma_total = sigma_0 + sigma_solidsolution + sigma_ppt;

    
        if (sigma_total > sigma_max) 
        {
            sigma_max = sigma_total;
            time_max = time;
        }

        time = time + dt;
        count += 1;
        
    }
    printf("time : %lf, vf: %lf, comp_matrix: %lf\n", time, vf, comp);
    printf("Sigma at final time: %lf MPa, Maximum sigma: %lf MPa, time_for maximum_sigma: %le s\n", sigma_total, sigma_max, time_max); 
    return 0;
}