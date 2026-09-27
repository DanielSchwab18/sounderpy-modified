#pragma once
#include <iostream>
#include <numbers>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <utility>
#include "Helper_Funcs.hpp"

class Interpolator {
public:
    Interpolator(std::string sounding_file_name);
    ~Interpolator();


private:
    double L_v = 2501000.0; //Latent heat of vaporization, using the value suggested in Peters et al 2023. Units are J/kg
    double R_v = 461.5; //Gas constant for water vapor. Units are J/(kg K)
    double T0 = 273.15; //T_0 for Clausius-Clapeyron. Units are K
    double e_s0 = 611.2; //e_{s0} for Clausius-Clapeyron. Units are Pa.
    double Lv_to_Rv_ratio = L_v / R_v; //L_v/R_v ratio. Units are K.
    double Rv_to_Lv_ratio = R_v / L_v; //R_v/L_v ratio. Units are K.
    double g = 9.81; //Constant gravity on earth, units are m/s^2
    double c_p = 1005.7;

    double R_d = 287.04; //Gas constant for dry air. Units are J/(kg K)
    double epsilon = R_d / R_v; //epsilon value. Unitless.
    double kappa = R_d / c_p;

    int interpolated_size;
    int size;

    double dry_adiabatic_lapse_rate = g / c_p;

    double max_interpolated_height = 0;

    double bunkers_u_peters = 0;
    double bunkers_v_peters = 0;
    double interpolated_bunkers_u = 0;
    double interpolated_bunkers_v = 0;

    double mixed_parcel_start_temp = 0;
    double mixed_parcel_start_dew = 0;
    double mixed_parcel_start_pressure = 0;

    double most_unstable_parcel_start_height = 0;
    double most_unstable_parcel_start_index = 0;

    double sb_lfc;
    int sb_lfc_index;
    double sb_el;
    int sb_el_index;
    double sb_cape;
    double sb_cin;
    double sb_ncape;
    double sb_ncape_peters;
    double sb_e_tilde;
    double sb_e_a;
    double sb_fractional_entrainment;
    double sb_ecape;
    double sb_ecin;
    double sb_e_lfc;
    double sb_e_el;
    double sb_3cape;
    double sb_3ecape;
        
    double ml_lfc;
    int ml_lfc_index;
    double ml_el;
    int ml_el_index;
    double ml_cape;
    double ml_cin;
    double ml_ncape;
    double ml_e_tilde;
    double ml_e_a;
    double ml_fractional_entrainment;
    double ml_ecape;
    double ml_ecin;
    double ml_e_lfc;
    double ml_e_el;
    double ml_3cape;
    double ml_3ecape;
    
    double mu_lfc;
    int mu_lfc_index;
    double mu_el;
    int mu_el_index;
    double mu_cape;
    double mu_cin;
    double mu_ncape;
    double mu_e_tilde;
    double mu_e_a;
    double mu_fractional_entrainment;
    double mu_ecape;
    double mu_ecin;
    double mu_e_lfc;
    double mu_e_el;
    double mu_3cape;
    double mu_3ecape;
    
    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> dewpoint;
    std::vector<double> relative_humidity;
    std::vector<double> wind_u;
    std::vector<double> wind_v;
    std::vector<double> height;
    std::vector<double> saturation_vapor_pressure;
    std::vector<double> vapor_pressure;
    std::vector<double> mixing_ratio;
    std::vector<double> specific_humidity;

    std::vector<double> interpolated_saturation_vapor_pressure;
    std::vector<double> interpolated_vapor_pressure;
    std::vector<double> interpolated_mixing_ratio;
    std::vector<double> lapse_rate; //lapse_rate[i] = lapse rate from height[i] to height[i+1]
    std::vector<double> interpolated_virtual_temperature;

    std::vector<double> interpolated_height;
    std::vector<double> interpolated_pressure;
    std::vector<double> interpolated_temperature;
    std::vector<double> interpolated_dewpoint;
    std::vector<double> interpolated_relative_humidity;
    std::vector<double> interpolated_wind_u;
    std::vector<double> interpolated_wind_v;
    std::vector<double> interpolated_specific_humidity;

    std::vector<double> storm_relative_wind_u_peters;
    std::vector<double> storm_relative_wind_v_peters;
    std::vector<double> storm_relative_wind_peters;
    
    std::vector<double> interpolated_storm_relative_wind_u;
    std::vector<double> interpolated_storm_relative_wind_v;
    std::vector<double> interpolated_storm_relative_wind;

    std::vector<double> interpolated_storm_relative_helicity;
    std::vector<double> storm_relative_helicity;

    std::vector<double> mu_parcel_buoyancy;
    std::vector<double> ml_parcel_buoyancy;
    std::vector<double> sb_parcel_buoyancy;
    
    std::vector<double> mu_entrainment_buoyancy;
    std::vector<double> ml_entrainment_buoyancy;
    std::vector<double> sb_entrainment_buoyancy;


    void calc_lapse_rates();

    void interpolate_temperature();
    void interpolate_pressure();
    void interpolate_relative_humidity();
    void interpolate_wind_u();
    void interpolate_wind_v();

    void calc_interpolated_saturation_vapor_pressures();
    void calc_interpolated_vapor_pressures();
    void calc_interpolated_mixing_ratios();
    void calc_interpolated_virtual_temperatures();

    void interpolate_dewpoint();

    void calc_right_moving_bunkers_motion_interpolated();
    void calc_right_moving_bunkers_motion_peters();
    void calc_storm_relative_winds();

    void calc_storm_relative_helicity();

    void calc_mixed_parcel_start_values();
    void calc_most_unstable_parcel_start_values();

    void calc_muparcel_path();
    void calc_mlparcel_path();
    void calc_sbparcel_path();

    double calc_saturation_mixing_ratio(double T, double p);
    double calc_moist_lapse_rate(double T, double p);

    void calc_sb_e_tilde_and_frac_entrainment();
    void calc_ml_e_tilde_and_frac_entrainment();
    void calc_mu_e_tilde_and_frac_entrainment();

    void calc_sb_entrainment_buoyancy();
    void calc_mu_entrainment_buoyancy();
    void calc_ml_entrainment_buoyancy();

    std::tuple<double, double, double, double, double> calc_cape_and_cin(std::vector<double> buoyancy_profile);

    double compute_mixing_ratio_saturation_linear_combo(double T, double p);
    double compute_mixing_ratio_saturation_all_ice(double T, double p);

    double dry_lapse_rate(double T, double qv, double T0, double qv0, double fracent);

    double calc_moist_lapse_rate_peters(double T, double qv, double qvv, double qvi, double p0, double T_0, double q0, double qt, double fracent);
    
    double calc_moist_lapse_rate_peters_no_e(double T, double qv, double qvv, double qvi, double p0, double T_0, double q0, double qt);

    double omega(double T);

    double domega(double T);

    void calc_sbncape();
    void calc_sbncape_peters();
    void calc_muncape();
    void calc_mlncape();
};

Sounding::Sounding(std::string sounding_file_name, bool interpolated_cape) {
    std::ifstream sounding_file(sounding_file_name);
    std::string curline;

    double p, t, td, rh, u, v, h;

    while (std::getline(sounding_file, curline)) {
        std::istringstream cur(curline);

        cur >> p >> t >> td >> rh >> u >> v >> h;

        pressure.push_back(p * 100);
        temperature.push_back(t + 273.15);
        dewpoint.push_back(td + 273.15);
        relative_humidity.push_back(rh);
        wind_u.push_back(u);
        wind_v.push_back(v);
        height.push_back(h);
    }

    sounding_file.close();
    
    for(int i = (int) height.size() - 1; i > 0; i--) {
        height[i] -= height[0];
    }
    height[0] = 0.0;

    calc_lapse_rates();

    int max_height = height[height.size() - 1];
    int interpolated_height_length = int(max_height / 100);
    for (int i = 0; i <= interpolated_height_length; i++) {
        interpolated_height.push_back(i * 100);
    }

    max_interpolated_height = interpolated_height[interpolated_height.size() - 1];
    interpolated_size = (int)interpolated_height.size();
    size = (int) height.size();

    interpolate_temperature();
    interpolate_relative_humidity();
    interpolate_wind_u();
    interpolate_wind_v();
    interpolate_pressure();

    for (int i = 0; i < interpolated_pressure.size() - 1; i++) {
        if (interpolated_pressure[i] < interpolated_pressure[i + 1]) {
            std::cout << "ERROR ERROR ERROR\n" << i << std::endl;
        }
    }

    calc_interpolated_saturation_vapor_pressures();
    calc_interpolated_vapor_pressures();
    interpolate_dewpoint();
    calc_interpolated_mixing_ratios();
    calc_interpolated_virtual_temperatures();

    
    calc_right_moving_bunkers_motion_interpolated();
    calc_right_moving_bunkers_motion_peters();
    calc_storm_relative_winds();

    calc_storm_relative_helicity();
    
    if(interpolated_cape) {
        pressure = interpolated_pressure;
        temperature = interpolated_temperature;
        height = interpolated_height;
        relative_humidity = interpolated_relative_humidity;
        saturation_vapor_pressure = interpolated_saturation_vapor_pressure;
        vapor_pressure = interpolated_vapor_pressure;
        mixing_ratio = interpolated_mixing_ratio;
        specific_humidity = interpolated_specific_humidity;
        size = interpolated_size;
    }
    
    calc_mixed_parcel_start_values();
    calc_most_unstable_parcel_start_values();

    calc_sbparcel_path();
    
    std::tuple<double, double, double, double, double> sb_out = calc_cape_and_cin(sb_parcel_buoyancy);
    sb_cape = std::get<0>(sb_out);
    sb_3cape = std::get<1>(sb_out);
    sb_cin = std::get<2>(sb_out);
    sb_lfc = std::get<3>(sb_out);
    sb_el = std::get<4>(sb_out);
    
    calc_sbncape();
    calc_sbncape_peters();
    calc_sb_e_tilde_and_frac_entrainment();
    calc_sb_entrainment_buoyancy();
    
    std::tuple<double, double, double, double, double> sb_e_out = calc_cape_and_cin(sb_entrainment_buoyancy);
    sb_ecape = std::get<0>(sb_e_out);
    sb_3ecape = std::get<1>(sb_e_out);
    sb_ecin = std::get<2>(sb_e_out);
    sb_e_lfc = std::get<3>(sb_e_out);
    sb_e_el = std::get<4>(sb_e_out);
    
    calc_mlparcel_path();
    
    std::tuple<double, double, double, double, double> ml_out = calc_cape_and_cin(ml_parcel_buoyancy);
    ml_cape = std::get<0>(ml_out);
    ml_3cape = std::get<1>(ml_out);
    ml_cin = std::get<2>(ml_out);
    ml_lfc = std::get<3>(ml_out);
    ml_el = std::get<4>(ml_out);
    
    calc_mlncape();
    calc_ml_e_tilde_and_frac_entrainment();
    calc_ml_entrainment_buoyancy();
    
    std::tuple<double, double, double, double, double> ml_e_out = calc_cape_and_cin(ml_entrainment_buoyancy);
    ml_ecape = std::get<0>(ml_e_out);
    ml_3ecape = std::get<1>(ml_e_out);
    ml_ecin = std::get<2>(ml_e_out);
    ml_e_lfc = std::get<3>(ml_e_out);
    ml_e_el = std::get<4>(ml_e_out);
    
    calc_muparcel_path();
    
    std::tuple<double, double, double, double, double> mu_out = calc_cape_and_cin(mu_parcel_buoyancy);
    mu_cape = std::get<0>(mu_out);
    mu_3cape = std::get<1>(mu_out);
    mu_cin = std::get<2>(mu_out);
    mu_lfc = std::get<3>(mu_out);
    mu_el = std::get<4>(mu_out);
    
    calc_muncape();
    calc_mu_e_tilde_and_frac_entrainment();
    calc_mu_entrainment_buoyancy();
    
    std::tuple<double, double, double, double, double> mu_e_out = calc_cape_and_cin(mu_entrainment_buoyancy);
    mu_ecape = std::get<0>(mu_e_out);
    mu_3ecape = std::get<1>(mu_e_out);
    mu_ecin = std::get<2>(mu_e_out);
    mu_e_lfc = std::get<3>(mu_e_out);
    mu_e_el = std::get<4>(mu_e_out);
}

Sounding::~Sounding() {}

void Sounding::calc_interpolated_saturation_vapor_pressures() {
    for (int i = 0; i < interpolated_temperature.size(); i++) {
        double t_diff = (1 / T0 - 1 / (interpolated_temperature[i]));
        double inner_exp = Lv_to_Rv_ratio * t_diff;
        double exponent = exp(inner_exp);
        interpolated_saturation_vapor_pressure.push_back(e_s0 * exponent);
    }
    
    for (int i = 0; i < temperature.size(); i++) {
        double t_diff = (1 / T0 - 1 / (temperature[i]));
        double inner_exp = Lv_to_Rv_ratio * t_diff;
        double exponent = exp(inner_exp);
        saturation_vapor_pressure.push_back(e_s0 * exponent);
    }
}

void Sounding::calc_interpolated_vapor_pressures() {
    for (int i = 0; i < interpolated_height.size(); i++) {
        interpolated_vapor_pressure.push_back(interpolated_saturation_vapor_pressure[i] * interpolated_relative_humidity[i] / 100);
    }
    
    for (int i = 0; i < height.size(); i++) {
        vapor_pressure.push_back(saturation_vapor_pressure[i] * relative_humidity[i] / 100);
    }
}

void Sounding::calc_interpolated_mixing_ratios() {
    for (int i = 0; i < interpolated_height.size(); i++) {
        //interpolated_mixing_ratio.push_back((epsilon * interpolated_vapor_pressure[i]) / (interpolated_pressure[i] - interpolated_vapor_pressure[i]));
        interpolated_mixing_ratio.push_back(compute_mixing_ratio_saturation_linear_combo(interpolated_dewpoint[i], interpolated_pressure[i]));
        interpolated_specific_humidity.push_back(interpolated_mixing_ratio[i] / (1 + interpolated_mixing_ratio[i]));
    }
    
    for (int i = 0; i < height.size(); i++) {
        mixing_ratio.push_back(compute_mixing_ratio_saturation_linear_combo(dewpoint[i], pressure[i]));
        specific_humidity.push_back(mixing_ratio[i] / (1 + mixing_ratio[i]));
    }
    
}

void Sounding::calc_lapse_rates() {
    for (int i = 0; i < height.size() - 1; i++) {
        double delta_T = temperature[i] - temperature[i + 1];
        double delta_z = height[i] - height[i + 1];
        lapse_rate.push_back(-1 * delta_T / delta_z);
    }
}

void Sounding::calc_interpolated_virtual_temperatures() {
    for (int i = 0; i < interpolated_height.size(); i++) {
        double q = interpolated_mixing_ratio[i] / (1 + interpolated_mixing_ratio[i]);
        double correction_factor = 1 + 0.608013937283 * q;
        interpolated_virtual_temperature.push_back(correction_factor * interpolated_temperature[i]);
    }
}

void Sounding::interpolate_temperature() {
    int cur_height_index = 1;
    if(std::fmod(height[height.size()-1], 100.0) == 0.0) {
        height[height.size()-1] += 1.0;
    }
    for (int i = 0; i < interpolated_height.size(); i++) {
        if (interpolated_height[i] < height[cur_height_index]) {
            double delta_z = height[cur_height_index] - interpolated_height[i];
            double delta_T = delta_z * lapse_rate[cur_height_index - 1];
            interpolated_temperature.push_back(temperature[cur_height_index] + delta_T);

        }
        else {
            cur_height_index++;
            i--;
        }
    }
}

void Sounding::interpolate_dewpoint() {
    for (int i = 0; i < interpolated_height.size(); i++) {
        double natural_log = std::log(interpolated_vapor_pressure[i] / e_s0);
        double subtract_part = Rv_to_Lv_ratio * natural_log;
        interpolated_dewpoint.push_back(1 / ((1 / T0) - subtract_part));
    }
}

void Sounding::interpolate_relative_humidity() {
    std::vector<double> change_in_rh;
    for (int i = 0; i < height.size() - 1; i++) {
        double delta_rh = relative_humidity[i] - relative_humidity[i + 1];
        double delta_z = height[i] - height[i + 1];
        change_in_rh.push_back(-1 * delta_rh / delta_z);
    }

    int cur_height_index = 1;
    for (int i = 0; i < interpolated_height.size(); i++) {
        if (interpolated_height[i] < height[cur_height_index]) {
            double delta_z = height[cur_height_index] - interpolated_height[i];
            double delta_rh = delta_z * change_in_rh[cur_height_index - 1];
            interpolated_relative_humidity.push_back(relative_humidity[cur_height_index] + delta_rh);
        }
        else {
            cur_height_index++;
            i--;
        }
    }
}

void Sounding::interpolate_wind_u() {
    std::vector<double> change_in_u;
    for (int i = 0; i < height.size() - 1; i++) {
        double delta_u = wind_u[i] - wind_u[i + 1];
        double delta_z = height[i] - height[i + 1];
        change_in_u.push_back(-1 * delta_u / delta_z);
    }

    int cur_height_index = 1;
    for (int i = 0; i < interpolated_height.size(); i++) {
        if (interpolated_height[i] < height[cur_height_index]) {
            double delta_z = height[cur_height_index] - interpolated_height[i];
            double delta_u = delta_z * change_in_u[cur_height_index - 1];
            interpolated_wind_u.push_back(wind_u[cur_height_index] + delta_u);
        }
        else {
            cur_height_index++;
            i--;
        }
    }
}

void Sounding::interpolate_wind_v() {
    std::vector<double> change_in_v;
    for (int i = 0; i < height.size() - 1; i++) {
        double delta_v = wind_v[i] - wind_v[i + 1];
        double delta_z = height[i] - height[i + 1];
        change_in_v.push_back(-1 * delta_v / delta_z);
    }

    int cur_height_index = 1;
    for (int i = 0; i < interpolated_height.size(); i++) {
        if (interpolated_height[i] < height[cur_height_index]) {
            double delta_z = height[cur_height_index] - interpolated_height[i];
            double delta_v = delta_z * change_in_v[cur_height_index - 1];
            interpolated_wind_v.push_back(wind_v[cur_height_index] + delta_v);
        }
        else {
            cur_height_index++;
            i--;
        }
    }
}

void Sounding::interpolate_pressure() {
    int cur_height_index = 1;
    for (int i = 0; i < interpolated_height.size(); i++) {
        if (interpolated_height[i] < height[cur_height_index]) {
            double exponent = g / (R_d * lapse_rate[cur_height_index - 1]);
            double fraction = (interpolated_temperature[i]) / (temperature[cur_height_index]);
            double factor = std::pow(fraction, exponent);
            interpolated_pressure.push_back(pressure[cur_height_index] * factor);

        }
        else {
            cur_height_index++;
            i--;
        }
    }
}

void Sounding::calc_right_moving_bunkers_motion_interpolated() {
    double u_sum = 0;
    double v_sum = 0;
    
    int i = 1;
    
    while(interpolated_height[i] < 6000) {
        u_sum += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        v_sum += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        i++;
    }
    if(interpolated_height[i] == 6000) {
        u_sum += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        v_sum += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
    } else {
        double height_diff = 6000 - interpolated_height[i-1];
        double height_prop = height_diff / (interpolated_height[i] - interpolated_height[i-1]);
        double interim_u = interpolated_wind_u[i-1] + (height_prop * (interpolated_wind_u[i] - interpolated_wind_u[i-1]));
        double interim_v = interpolated_wind_v[i-1] + (height_prop * (interpolated_wind_v[i] - interpolated_wind_v[i-1]));
        u_sum += (interim_u + interpolated_wind_u[i-1]) / 2 * height_diff;
        v_sum += (interim_v + interpolated_wind_v[i-1]) / 2 * height_diff;
    }

    double mean_u = u_sum / 6000;
    double mean_v = v_sum / 6000;

    double bottom_500_mean_u = 0;
    double bottom_500_mean_v = 0;
    double top_500_mean_u = 0;
    double top_500_mean_v = 0;
    
    i=1;
    while(interpolated_height[i] < 500) {
        bottom_500_mean_u += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        bottom_500_mean_v += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        i++;
    }
    if(interpolated_height[i] == 500) {
        bottom_500_mean_u += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        bottom_500_mean_v += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
    } else {
        double height_diff = 500 - interpolated_height[i-1];
        double height_prop = height_diff / (interpolated_height[i] - interpolated_height[i-1]);
        double interim_u = interpolated_wind_u[i-1] + (height_prop * (interpolated_wind_u[i] - interpolated_wind_u[i-1]));
        double interim_v = interpolated_wind_v[i-1] + (height_prop * (interpolated_wind_v[i] - interpolated_wind_v[i-1]));
        bottom_500_mean_u += (interim_u + interpolated_wind_u[i-1]) / 2 * height_diff;
        bottom_500_mean_v += (interim_v + interpolated_wind_v[i-1]) / 2 * height_diff;
    }
    
    while(interpolated_height[i] < 5500) {
        i++;
    }
    
    if(interpolated_height[i] != 5500) {
        double height_diff = interpolated_height[i] - 5500;
        double height_prop = height_diff / (interpolated_height[i] - interpolated_height[i-1]);
        double interim_u = interpolated_wind_u[i] - (height_prop * (interpolated_wind_u[i] - interpolated_wind_u[i-1]));
        double interim_v = interpolated_wind_v[i] - (height_prop * (interpolated_wind_v[i] - interpolated_wind_v[i-1]));
        top_500_mean_u += (interim_u + interpolated_wind_u[i]) / 2 * height_diff;
        top_500_mean_v += (interim_v + interpolated_wind_v[i]) / 2 * height_diff;
    }
    i++;
    while(interpolated_height[i] < 6000) {
        top_500_mean_u += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        top_500_mean_v += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        i++;
    }
    if(interpolated_height[i] == 6000) {
        top_500_mean_u += (interpolated_wind_u[i] + interpolated_wind_u[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
        top_500_mean_v += (interpolated_wind_v[i] + interpolated_wind_v[i-1]) / 2 * (interpolated_height[i] - interpolated_height[i-1]);
    } else {
        double height_diff = 6000 - interpolated_height[i-1];
        double height_prop = height_diff / (interpolated_height[i] - interpolated_height[i-1]);
        double interim_u = interpolated_wind_u[i-1] + (height_prop * (interpolated_wind_u[i] - interpolated_wind_u[i-1]));
        double interim_v = interpolated_wind_v[i-1] + (height_prop * (interpolated_wind_v[i] - interpolated_wind_v[i-1]));
        top_500_mean_u += (interim_u + interpolated_wind_u[i-1]) / 2 * height_diff;
        top_500_mean_v += (interim_v + interpolated_wind_v[i-1]) / 2 * height_diff;
    }

    bottom_500_mean_u /= 500.0;
    bottom_500_mean_v /= 500.0;
    top_500_mean_u /= 500.0;
    top_500_mean_v /= 500.0;

    double bunkers_rm_u_change = top_500_mean_v - bottom_500_mean_v;
    double bunkers_rm_v_change = bottom_500_mean_u - top_500_mean_u;

    double length = std::pow(bunkers_rm_u_change * bunkers_rm_u_change + bunkers_rm_v_change * bunkers_rm_v_change, 0.5);

    bunkers_rm_u_change /= length;
    bunkers_rm_v_change /= length;
    
    interpolated_bunkers_u = mean_u + bunkers_rm_u_change * 7.5;
    interpolated_bunkers_v = mean_v + bunkers_rm_v_change * 7.5;

}

void Sounding::calc_right_moving_bunkers_motion_peters() {
    double u_sum = 0;
    double v_sum = 0;
    double bottom_500_mean_u = 0;
    double bottom_500_mean_v = 0;
    double top_500_mean_u = 0;
    double top_500_mean_v = 0;
    
    int i = 0;
    int bottomcounts = 0;
    int topcounts = 0;
    
    while(height[i] <= 6000) {
        u_sum += wind_u[i] * (height[i+1] - height[i]);
        v_sum += wind_v[i] * (height[i+1] - height[i]);
        if(height[i] <= 500) {
            bottom_500_mean_u += wind_u[i];
            bottom_500_mean_v += wind_v[i];
            bottomcounts++;
        } else if (height[i] >= 5500) {
            top_500_mean_u += wind_u[i];
            top_500_mean_v += wind_v[i];
            topcounts++;
        }
        i++;
    }


    double mean_u = u_sum / height[i];
    double mean_v = v_sum / height[i];

    bottom_500_mean_u /= bottomcounts;
    bottom_500_mean_v /= bottomcounts;
    top_500_mean_u /= topcounts;
    top_500_mean_v /= topcounts;

    double bunkers_rm_u_change = top_500_mean_v - bottom_500_mean_v;
    double bunkers_rm_v_change = bottom_500_mean_u - top_500_mean_u;
    
    double length = std::pow(bunkers_rm_u_change * bunkers_rm_u_change + bunkers_rm_v_change * bunkers_rm_v_change, 0.5);

    bunkers_rm_u_change /= length;
    bunkers_rm_v_change /= length;
    
    std::vector<double> dudz;
    std::vector<double> dvdz;
    for (int i = 1; i < height.size() - 1; i++) {
        dudz.push_back((wind_u[i+1] - wind_u[i-1]) / (height[i+1] - height[i-1]));
        dvdz.push_back((wind_v[i+1] - wind_v[i-1]) / (height[i+1] - height[i-1]));
    }
    dudz.insert(dudz.begin(), 2 * dudz[0] - dudz[1]);
    dvdz.insert(dvdz.begin(), 2 * dvdz[0] - dvdz[1]);
    i = 0;
    double sum = 0;
    while(height[i] <= 1000) {
        sum += -1 * (wind_u[i] - mean_u) * dvdz[i] + (wind_v[i] - mean_v) * dudz[i];
        i++;
    }
    sum *= 1000;
    double SRH_mean = sum/i;
    double propfac = 0.0;
    if(SRH_mean < -150 || SRH_mean > 150) {
        propfac = 1.0;
    } else if (SRH_mean == 0) {
        propfac = 0.0;
    } else if (SRH_mean < 0) {
        propfac = -1 * SRH_mean / 150.0;
    } else {
        propfac = SRH_mean / 150.0;
    }
    
    bunkers_u_peters = mean_u + propfac * (bunkers_rm_u_change * 7.5);
    bunkers_v_peters = mean_v + propfac * (bunkers_rm_v_change * 7.5);
}


void Sounding::calc_storm_relative_winds() {
    for (int i = 0; i < interpolated_wind_u.size(); i++) {
        interpolated_storm_relative_wind_u.push_back(interpolated_wind_u[i] - interpolated_bunkers_u);
        interpolated_storm_relative_wind_v.push_back(interpolated_wind_v[i] - interpolated_bunkers_v);
        interpolated_storm_relative_wind.push_back(std::pow(std::pow(interpolated_storm_relative_wind_u[i], 2.0) + std::pow(interpolated_storm_relative_wind_v[i], 2.0), 0.5));
    }
    
    for (int i = 0; i < wind_u.size(); i++) {
        storm_relative_wind_u_peters.push_back(wind_u[i] - bunkers_u_peters);
        storm_relative_wind_v_peters.push_back(wind_v[i] - bunkers_v_peters);
        storm_relative_wind_peters.push_back(std::pow(std::pow(storm_relative_wind_u_peters[i], 2.0) + std::pow(storm_relative_wind_v_peters[i], 2.0), 0.5));
    }
}

void Sounding::calc_storm_relative_helicity() {
    for (int i = 0; i < interpolated_height.size() - 1; i++) {
        double first_term = (interpolated_wind_u[i + 1] - interpolated_bunkers_u) * (interpolated_wind_v[i] - interpolated_bunkers_v);
        double second_term = (interpolated_wind_u[i] - interpolated_bunkers_u) * (interpolated_wind_v[i + 1] - interpolated_bunkers_v);;
        
        double srh = first_term - second_term;
        interpolated_storm_relative_helicity.push_back(srh);
    }
    
    for (int i = 0; i < height.size() - 1; i++) {
        double first_term = (wind_u[i + 1] - bunkers_u_peters) * (wind_v[i] - bunkers_v_peters);
        double second_term = (wind_u[i] - bunkers_u_peters) * (wind_v[i + 1] - bunkers_v_peters);;
        
        double srh = first_term - second_term;
        storm_relative_helicity.push_back(srh);
    }
    
}

void Sounding::calc_mixed_parcel_start_values() {
    double ground_pressure = interpolated_pressure[0];
    double pressure_bound = ground_pressure - 10000.0;

    int i = 0;
    double pressure_sum = 0;
    double temperature_sum = 0;
    double dewpoint_sum = 0;
    while (interpolated_pressure[i] >= pressure_bound) {
        pressure_sum += interpolated_pressure[i];
        temperature_sum += interpolated_temperature[i];
        dewpoint_sum += interpolated_dewpoint[i];
        i++;
    }

    mixed_parcel_start_temp = temperature_sum / i;
    mixed_parcel_start_dew = dewpoint_sum / i;
    mixed_parcel_start_pressure = pressure_sum / i;
}

void Sounding::calc_most_unstable_parcel_start_values() {
    double max_moist_static_energy = 0;
    for (int i = 0; pressure[i] >= pressure[0]-30000; i++) {
        double MSE = c_p * temperature[i] + L_v * specific_humidity[i] + g * height[i];
        if(MSE > max_moist_static_energy) {
            max_moist_static_energy = MSE;
            most_unstable_parcel_start_height = height[i];
            most_unstable_parcel_start_index = i;
        }
    }
}

double Sounding::calc_saturation_mixing_ratio(double T, double p) {
    double cpl = 4190.0;
    double cpv = 1870.0;
    double term_one = (cpv - cpl) / R_v;
    double term_two = (L_v - T0 * (cpv - cpl)) / R_v;
    double esl = std::exp((T - T0) * term_two / (T * T0)) * e_s0 * std::pow((T / T0), term_one);

    return epsilon * esl / (p - esl);
}

double Sounding::calc_moist_lapse_rate(double T, double p) {
    double w_s = calc_saturation_mixing_ratio(T, p);
    double frac_one = (L_v * w_s) / (R_d * T);
    double frac_two = (L_v * L_v * w_s) / (R_v * c_p * T * T);
    double part_one = (1.0 + frac_one) / (1.0 + frac_two);
    return (part_one * dry_adiabatic_lapse_rate);
}

void Sounding::calc_sbparcel_path() {
    //This is functionally identical to the code from Peters 2023.
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    T_lif.push_back(temperature[0]);
    Q_v_lif.push_back(specific_humidity[0]);
    Q_t_lif.push_back(specific_humidity[0]);
    
    //Helper_Funcs::printVector(specific_humidity);
    
    double q_sat_prev = 0;
    double i = 0;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if (Q_v_lif[i-1]<q_sat) {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0));
            //std::cout << "i: " << i << " T_lif[i]: " << T_lif[i] << std::endl;
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            
            if (Q_v_lif[i] >= q_sat) {
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);



                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0);
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + pressure[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters_no_e(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters_no_e(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1]));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (0 * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));

            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    
    /*
    std::cout << "T_lif" << std::endl;
    Helper_Funcs::printVector(T_lif);
    std::cout << "Q_t_lif" << std::endl;
    Helper_Funcs::printVector(Q_t_lif);
    std::cout << "Q_v_lif" << std::endl;
    Helper_Funcs::printVector(Q_v_lif);*/
    
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v/R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v/R_d - 1) * specific_humidity[i]));
    }
    /*
    std::cout << "T_rho_lif" << std::endl;
    Helper_Funcs::printVector(T_rho_lif);
    std::cout << "T_0_lif" << std::endl;
    Helper_Funcs::printVector(T_0_lif);*/
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        sb_parcel_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
    
    //Helper_Funcs::printVector(T_rho_lif);
    
}

//returns {CAPE, 3CAPE, CIN, LFC, EL}
std::tuple<double, double, double, double, double> Sounding::calc_cape_and_cin(std::vector<double> buoyancy_profile) {
    int i = (int)buoyancy_profile.size() - 1;
    double CAPE = 0;
    double CAPE3k = 0;
    double CIN = 0;
    double LFC = 0;
    double EL = 0;
    while (buoyancy_profile[i] < 0) {
        if(buoyancy_profile[i-1] >= 0) {
            double buoyancy_change = buoyancy_profile[i-1] - buoyancy_profile[i];
            double height_ratio = buoyancy_profile[i-1] / buoyancy_change;
            EL = height_ratio * (height[i] - height[i-1]) + height[i-1];
            if(EL <= 3000) {
                CAPE3k += (buoyancy_profile[i-1]) / 2 * (EL - height[i-1]);
            }
            CAPE += (buoyancy_profile[i-1]) / 2 * (EL - height[i-1]);
        }
        i--;
        if(i == 0) {
            return {CAPE, CAPE3k, CIN, LFC, EL};
        }
    }
    bool keepGoing = true;
    while (keepGoing) {
        if(buoyancy_profile[i-1] <= 0) {
            if(buoyancy_profile[i] < 0) {
                i--;
                continue;
            }
            double buoyancy_change = buoyancy_profile[i] - buoyancy_profile[i-1];
            double height_ratio = -1 * buoyancy_profile[i-1] / buoyancy_change;
            LFC = height_ratio * (height[i] - height[i-1]) + height[i-1];
            if(height[i] <= 3000) {
                CAPE3k += (buoyancy_profile[i]) / 2 * (height[i] - LFC);
            }
            CAPE += (buoyancy_profile[i]) / 2 * (height[i] - LFC);
            CIN += (buoyancy_profile[i-1]) / 2 * (LFC - height[i-1]);
            
            keepGoing = false;
            for(int j = i-1; j > 1; j--) {
                if(buoyancy_profile[j] > 0) {
                    keepGoing = true;
                }
            }
        } else {
            if(height[i] <= 3000) {
                CAPE3k += (buoyancy_profile[i] + buoyancy_profile[i-1]) / 2 * (height[i] - height[i-1]);
            }
            CAPE += (buoyancy_profile[i] + buoyancy_profile[i-1]) / 2 * (height[i] - height[i-1]);
        }
        
        i--;
        if(i == 0) {
            break;
        }
    }
    
    while (i > 0) {
        if(buoyancy_profile[i] > 0) {
            i--;
            continue;
        }
        if(i == 1) {
            CIN += height[i] * buoyancy_profile[i] / 2;
        } else {
            CIN += (buoyancy_profile[i] + buoyancy_profile[i-1]) / 2 * (height[i] - height[i-1]);
        }
        i--;
    }
    
    return {CAPE, CAPE3k, CIN, LFC, EL};
}

void Sounding::calc_sbncape() {
    std::vector<double> moist_static_energy;
    std::vector<double> saturated_moist_static_energy;
    for (int i = 0; i < size; i++) {
        moist_static_energy.push_back(c_p * temperature[i] + L_v * specific_humidity[i] + g * height[i]);

        double saturation_mixing_ratio = calc_saturation_mixing_ratio(temperature[i], pressure[i]);
        double qsat = (1 - saturation_mixing_ratio) * saturation_mixing_ratio;
        //std::cout << "rsat: " << saturation_mixing_ratio << ", qsat: " << qsat << std::endl;
        saturated_moist_static_energy.push_back(c_p * temperature[i] + L_v * qsat + g * height[i]);
    }
    
    /*std::cout << "q0:" << std::endl;
    Helper_Funcs::printVector(specific_humidity);
    std::cout << "MSE0:" << std::endl;
    Helper_Funcs::printVector(moist_static_energy);*/

    std::vector<double> moist_static_energy_bar;
    moist_static_energy_bar.push_back(moist_static_energy[0]);

    std::vector<double> MSE_trap_sum;
    std::vector<double> height_trap_sum;
    std::vector<double> MSE_times_height;

    for (int i = 0; i < size - 1; i++) {
        MSE_trap_sum.push_back(moist_static_energy[i] + moist_static_energy[i + 1]);
        height_trap_sum.push_back(height[i + 1] - height[i]);
    }

    for (int i = 0; i < MSE_trap_sum.size(); i++) {
        MSE_times_height.push_back(MSE_trap_sum[i] * height_trap_sum[i]);
    }

    for (int i = 1; i < size; i++) {
        double numerator = Helper_Funcs::sumVector(MSE_times_height, 0, i);
        double denominator = 2 * (height[i] - height[0]);
        moist_static_energy_bar.push_back(numerator / denominator);
    }

    
    
    std::vector<double> int_arg;

    for (int i = 0; i < size; i++) {
        double term_one = g / (c_p * temperature[i]);
        double term_two = moist_static_energy_bar[i] - saturated_moist_static_energy[i];
        int_arg.push_back(-1.0 * term_one * term_two);
    }
    

    double sum = 0.0;
    int j = 0;
    for (int i = 0; height[i+1] < sb_el; i++) {
        if(height[i] < sb_lfc) {
            if(height[i+1] < sb_lfc) {
                continue;
            }
            double height_diff = (height[i+1] - sb_lfc);
            double height_prop = height_diff / (height[i+1] - height[i]);
            double int_arg_diff = int_arg[i+1] - int_arg[i];
            double interim_int_arg = (int_arg_diff * height_prop) + int_arg[i];
            double term_one = 0.5 * (interim_int_arg + int_arg[i + 1]);
            double term_two = height[i + 1] - sb_lfc;
            sum += term_one * term_two;
            continue;
            
        }
        double term_one = 0.5 * (int_arg[i] + int_arg[i + 1]);
        double term_two = height[i + 1] - height[i];
        sum += term_one * term_two;
        j = i;
    }
    j++;
    double height_diff = (height[j+1] - sb_el);
    double height_prop = height_diff / (height[j+1] - height[j]);
    double int_arg_diff = int_arg[j+1] - int_arg[j];
    double interim_int_arg = (int_arg_diff * height_prop) + int_arg[j];
    double term_one = 0.5 * (int_arg[j] + interim_int_arg);
    double term_two = sb_el - height[j];
    sum += term_one * term_two;
    
    
    
    /*std::cout << "int_arg" << std::endl;
    Helper_Funcs::printVector(int_arg);
    std::cout << "MSE0_star" << std::endl;
    Helper_Funcs::printVector(saturated_moist_static_energy);
    std::cout << "MSE0bar" << std::endl;
    Helper_Funcs::printVector(moist_static_energy_bar);
    std::cout << "T0" << std::endl;
    Helper_Funcs::printVector(temperature);
    std::cout << "lfc: " << sb_lfc_index << ", el: " << sb_el_index << std::endl;
    std::cout << "sum: " << sum << std::endl;
    */
    
    if(sum < 0) {
        sum = 0;
    }
    
    sb_ncape = sum;
}

void Sounding::calc_sbncape_peters() {
    std::vector<double> moist_static_energy;
    std::vector<double> saturated_moist_static_energy;
    for (int i = 0; i < size; i++) {
        moist_static_energy.push_back(c_p * temperature[i] + L_v * specific_humidity[i] + g * height[i]);

        double saturation_mixing_ratio = calc_saturation_mixing_ratio(temperature[i], pressure[i]);
        double qsat = (1 - saturation_mixing_ratio) * saturation_mixing_ratio;
        //std::cout << "rsat: " << saturation_mixing_ratio << ", qsat: " << qsat << std::endl;
        saturated_moist_static_energy.push_back(c_p * temperature[i] + L_v * qsat + g * height[i]);
    }
    
    /*std::cout << "q0:" << std::endl;
    Helper_Funcs::printVector(specific_humidity);
    std::cout << "MSE0:" << std::endl;
    Helper_Funcs::printVector(moist_static_energy);*/

    std::vector<double> moist_static_energy_bar;
    moist_static_energy_bar.push_back(moist_static_energy[0]);

    std::vector<double> MSE_trap_sum;
    std::vector<double> height_trap_sum;
    std::vector<double> MSE_times_height;

    for (int i = 0; i < size - 1; i++) {
        MSE_trap_sum.push_back(moist_static_energy[i] + moist_static_energy[i + 1]);
        height_trap_sum.push_back(height[i + 1] - height[i]);
    }

    for (int i = 0; i < MSE_trap_sum.size(); i++) {
        MSE_times_height.push_back(MSE_trap_sum[i] * height_trap_sum[i]);
    }

    for (int i = 1; i < size; i++) {
        double numerator = Helper_Funcs::sumVector(MSE_times_height, 0, i);
        double denominator = 2 * (height[i] - height[0]);
        moist_static_energy_bar.push_back(numerator / denominator);
    }

    
    
    std::vector<double> int_arg; //I have honestly no idea what this stands for, but it's the abbreviation Peters used.

    for (int i = 0; i < size; i++) {
        double term_one = g / (c_p * temperature[i]);
        double term_two = moist_static_energy_bar[i] - saturated_moist_static_energy[i];
        int_arg.push_back(-1.0 * term_one * term_two);
    }
    
    
    
    double sum = 0.0;
    for (int i = 1; height[i+1] <= sb_el; i++) {
        if(height[i] < sb_lfc) {
            continue;
        }
        sum += (int_arg[i] + int_arg[i + 1]) * (height[i+1] - height[i]) / 2;
    }
        
    
    /*
    std::cout << "int_arg" << std::endl;
    Helper_Funcs::printVector(int_arg);
    std::cout << "MSE0_star" << std::endl;
    Helper_Funcs::printVector(saturated_moist_static_energy);
    std::cout << "MSE0bar" << std::endl;
    Helper_Funcs::printVector(moist_static_energy_bar);
    std::cout << "T0" << std::endl;
    Helper_Funcs::printVector(temperature);
    std::cout << "lfc: " << sb_lfc_index << ", el: " << sb_el_index << std::endl;
    std::cout << "sum: " << sum << std::endl;
    */
    
    if(sum < 0) {
        sum = 0;
    }
    
    sb_ncape_peters = sum;
}


void Sounding::calc_sb_e_tilde_and_frac_entrainment() {
    double v_sr = get_srw(0, 1000);
    double H = sb_el;
    double sigma = 1.1;
    double alpha = 0.8;
    double prandtl_number = 1.0 / 3.0;
    double von_karman_constant = 0.18;
    double psi = von_karman_constant * (alpha * alpha) * (std::numbers::pi * std::numbers::pi) * 120.0 / (4.0 * prandtl_number * (sigma * sigma) * H);
    double vsr_tilde = v_sr / (std::pow(2.0 * sb_cape, 0.5));
    double n_tilde = sb_ncape / sb_cape;
    
    double vsr_tilde_squared = vsr_tilde * vsr_tilde;

    double numerator = -1 - psi - (psi / vsr_tilde_squared) * n_tilde;
    double root_term_one = 1 + psi + (psi / vsr_tilde_squared) * n_tilde;
    double root_term_two = 4 * psi / vsr_tilde_squared * (1 - psi * n_tilde);
    double square_root = std::pow(root_term_one * root_term_one + root_term_two, 0.5);

    numerator += square_root;

    double fraction = numerator / (2 * psi / vsr_tilde_squared);
    sb_e_tilde = fraction;

    double e_tilde = fraction;

    sb_fractional_entrainment = 2 * ((1 - e_tilde) / (e_tilde + n_tilde)) / (sb_el);
    
    sb_e_a = vsr_tilde_squared + (-1 - psi - (psi / vsr_tilde_squared) * n_tilde + std::pow(std::pow(1 + psi + (psi / vsr_tilde_squared) * n_tilde, 2) + 4 * (psi / vsr_tilde_squared) * (1 - psi * n_tilde), 0.5)) / (2 * psi / (vsr_tilde_squared));
    
    //std::cout << e_tilde << ", " << sb_e_tilde << ", " << n_tilde << ", " << sb_fractional_entrainment << std::endl;
}

void Sounding::calc_sb_entrainment_buoyancy() {
    //This is functionally identical to the code from Peters 2023.
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    T_lif.push_back(temperature[0]);
    Q_v_lif.push_back(specific_humidity[0]);
    Q_t_lif.push_back(specific_humidity[0]);

    double q_sat_prev = 0;
    double i = 0;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if ((Q_v_lif[i - 1] < q_sat) && (pressure[i] > 50000)){
            //std::cout << "i: " << i << std::endl;
            //std::cout << "qsat: " << q_sat << std::endl;
            //std::cout << "Qv_lif[i-1]: " << Q_v_lif[i-1] << std::endl;
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], sb_fractional_entrainment));
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * sb_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            //std::cout << Q_v_lif[i] << std::endl;
            //std::cout << q_sat << std::endl;
            //std::cout << q_sat_prev << std::endl;
            
            if (Q_v_lif[i] >= q_sat) {
                
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                
                if(satrat < 0.05) {
                    q_sat_prev = q_sat;
                    continue;
                }
                
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);
                //std::cout << "satrat: " << satrat << std::endl;


                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], sb_fractional_entrainment);
                //std::cout << "T_halfstep: " << T_halfstep << std::endl;
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * sb_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + specific_humidity[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep, sb_fractional_entrainment);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * sb_fractional_entrainment * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1], sb_fractional_entrainment));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (sb_fractional_entrainment * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));
            
            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }
    
    

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v / R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v / R_d - 1) * specific_humidity[i]));
    }
    
    /*
    std::cout << "T_lif: " << std::endl;
    Helper_Funcs::printVector(T_lif);
    std::cout << "Q_v_lif: " << std::endl;
    Helper_Funcs::printVector(Q_v_lif);
    std::cout << "Q_t_lif: " << std::endl;
    Helper_Funcs::printVector(Q_t_lif);
    std::cout << "T_rho_lif: " << std::endl;
    Helper_Funcs::printVector(T_rho_lif);
    std::cout << "T_0_lif: " << std::endl;
    Helper_Funcs::printVector(T_0_lif);*/
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        sb_entrainment_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
}

double Sounding::compute_mixing_ratio_saturation_linear_combo(double T, double p) {
    //This is identical to option 1 (linear combo) of Peters' 'compute_rsat' function.

    double cpv = 1870.0;
    double cpl = 4190.0;
    double cpi = 2106.0;
    double xls = 2834000.0;

    double omeg = omega(T);

    double term1 = (cpv - cpl) / R_v;
    double term2 = (L_v - T0 * (cpv - cpl)) / R_v;
    double esl_l = std::exp((T - T0) * term2 / (T * T0)) * e_s0 * std::pow((T / T0), term1);
    double qsat_l = epsilon * esl_l / (p - esl_l);
    term1 = (cpv - cpi) / R_v;
    term2 = (xls - T0 * (cpv - cpi)) / R_v;
    double esl_i = std::exp((T - T0) * term2 / (T * T0)) * e_s0 * std::pow((T / T0), term1);
    double qsat_i = epsilon * esl_i / (p - esl_i);
    return  (1.0 - omeg) * qsat_l + (omeg)*qsat_i;
}

double Sounding::dry_lapse_rate(double T, double qv, double T0, double qv0, double fracent) {
    //Functionally identical to Peters' drylift function.
    double cpv = 1870.0;

    double cpmv = (1.0 - qv) * c_p + qv * cpv;
    double B = g * ((T - T0) / T0 + (R_v / R_d - 1) * (qv - qv0));
    double eps = -1.0 * fracent * (T - T0);
    
    //std::cout << T << "," << T0 << "," << qv << "," << qv0 << std::endl;
    
    return -1 * ((g + B) / cpmv + eps);
}

double Sounding::compute_mixing_ratio_saturation_all_ice(double T, double p) {
    double cpv = 1870.0;
    double cpi = 2106.0;
    double xls = 2834000.0;

    double term1 = (cpv - cpi) / R_v;
    double term2 = (xls - T0 * (cpv - cpi)) / R_v;
    double esl = std::exp((T - T0) * term2 / (T * T0)) * e_s0 * std::pow((T / T0), term1);
    if (p * 0.5 < esl) {
        esl = p * 0.5;
    }
    return epsilon * esl / (p - esl);
}

double Sounding::calc_moist_lapse_rate_peters(double T, double qv, double qvv, double qvi, double p0, double T_0, double q0, double qt, double fracent) {
    double xls = 2834000.0;
    double cpv = 1870.0;
    double cpl = 4190.0;
    double cpi = 2106.0;
    
    double zeroB = 1;
    bool noT = false;
    double pseudofac = 1 - qv;
    if(qt == qv) {
        pseudofac = 1;
    }
    
    if (qt < 0) {
        qt = -1 * qt;
        zeroB = 0;
    }
    if (qv < 0) {
        qv = -1 * qv;
        noT = true;
    }
    
    double OMEGA = omega(T);
    double dOMEGA = domega(T);


    double cpm = (1 - qt) * c_p + qv * cpv + (1 - OMEGA) * (qt - qv) * cpl + OMEGA * (qt - qv) * cpi;
    double Lv = L_v + (T - T0) * (cpv - cpl);
    double Li = (xls - L_v) + (T - T0) * (cpl - cpi);
    double T_rho = T * (1.0 - qt + qv / epsilon);

    double Rm0 = (1 - q0) * R_d + q0 * R_v;
    double T_rho0 = T_0 * (1 - q0 + q0 / epsilon);
    double B = zeroB * g * (T_rho - T_rho0) / (T_rho0);

    double eps_T = -1.0 * fracent * (T - T_0);
    double eps_qv = -1.0 * fracent * (qv - q0);
    double eps_qt = -1.0 * fracent * (qt - q0);
   

    double Qvsl = qvv / (epsilon - epsilon * qt + qv);
    double Qvsi = qvi / (epsilon - epsilon * qt + qv);
    double Q_M = (1 - OMEGA) * qvv / (1 - Qvsl) + OMEGA * qvi / (1 - Qvsi);
    double L_M = Lv * (1 - OMEGA) * qvv / (1 - Qvsl) + (Lv + Li) * OMEGA * qvi / (1 - Qvsi);


    double term1 = -1 * B;

    double term2 = 0;
    if (noT) {
        term2 = -1 * pseudofac * Q_M * (Lv + Li * OMEGA) * g / (R_d * T_rho);
    } else {
        term2 = -1 * pseudofac * Q_M * (Lv + Li * OMEGA) * g / (Rm0 * T_0);
    }
    double term3 = -1 * g;
    double term4 = (cpm - Li * (qt - qv) * dOMEGA) * eps_T;
    double term5 = (Lv + Li * OMEGA) * (eps_qv + (qv / (1 - qt)) * eps_qt);

    double term6 = cpm;
    double term7 = -1 * Li * (qt - qv) * dOMEGA;
    double term8 = (Lv + Li * OMEGA) * (-1 * dOMEGA * (qvv - qvi) + (1 / (R_v * (T*T))) * (L_M));
    double gamma_m = (term1 + term2 + term3 + term4 + term5) / (term6 + term7 + term8);
    return gamma_m;
    
}

double Sounding::calc_moist_lapse_rate_peters_no_e(double T, double qv, double qvv, double qvi, double p0, double T_0, double q0, double qt ) {
    double xls = 2834000.0;
    double cpv = 1870.0;
    double cpl = 4190.0;
    double cpi = 2106.0;
    
    if (qt < 0) {
        qt = 0;
    }
    if (qv < 0) {
        qv = 0;
    }
    
    double OMEGA = omega(T);
    double dOMEGA = domega(T);


    double cpm = (1 - qt) * c_p + qv * cpv + (1 - OMEGA) * (qt - qv) * cpl + OMEGA * (qt - qv) * cpi;
    double Lv = L_v + (T - T0) * (cpv - cpl);
    double Li = (xls - L_v) + (T - T0) * (cpl - cpi);
    double T_rho = T * (1.0 - qt + qv / epsilon);




    double Rm0 = (1 - q0) * R_d + q0 * R_v;
    double T_rho0 = T_0 * (1 - q0 + q0 / epsilon);
    double B = g * (T_rho - T_rho0) / (T_rho0);

    double eps_T = -1.0 * 0 * (T - T_0);
    double eps_qv = -1.0 * 0 * (qv - q0);
    double eps_qt = -1.0 * 0 * (qt - q0);
   

    double Qvsl = qvv / (epsilon - epsilon * qt + qv);
    double Qvsi = qvi / (epsilon - epsilon * qt + qv);
    double Q_M = (1 - OMEGA) * qvv / (1 - Qvsl) + OMEGA * qvi / (1 - Qvsi);
    double L_M = Lv * (1 - OMEGA) * qvv / (1 - Qvsl) + (Lv + Li) * OMEGA * qvi / (1 - Qvsi);


    double term1 = -B;

    double term2 = -Q_M * (Lv + Li * OMEGA) * g / (Rm0 * T_0);

    double term3 = -g;
    double term4 = (cpm - Li * (qt - qv) * dOMEGA) * eps_T;
    double term5 = (Lv + Li * OMEGA) * (eps_qv + (qv / (1 - qt)) * eps_qt);

    double term6 = cpm;
    double term7 = -Li * (qt - qv) * dOMEGA;
    double term8 = (Lv + Li * OMEGA) * (-dOMEGA * (qvv - qvi) + (1 / (R_v * (T*T))) * (L_M));
    double gamma_m = (term1 + term2 + term3 + term4 + term5) / (term6 + term7 + term8);
    return gamma_m;
    
}



double Sounding::omega(double T) {
    if(T >= 273.15) {
        return 0.0;
    }
    if(T <= 253.15) {
        return 1;
    }
    
    return ((273.15 - T) / 20);
}

double Sounding::domega(double T) {
    if(T > 253.15 && T <= 273.15) {
        return -0.05;
    }
    return 0;
}

void Sounding::calc_mlparcel_path() {
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    T_lif.push_back(mixed_parcel_start_temp);
    double starting_mixing_ratio = calc_saturation_mixing_ratio(mixed_parcel_start_dew, mixed_parcel_start_pressure);
    double starting_specific_humidity = starting_mixing_ratio / (1 + starting_mixing_ratio);
    Q_v_lif.push_back(starting_specific_humidity);
    Q_t_lif.push_back(starting_specific_humidity);
        
    double q_sat_prev = 0;
    double i = 0;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if (Q_v_lif[i-1]<q_sat) {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0));
            //std::cout << "i: " << i << " T_lif[i]: " << T_lif[i] << std::endl;
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            
            if (Q_v_lif[i] >= q_sat) {
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);



                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0);
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + pressure[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters_no_e(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters_no_e(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1]));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (0 * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));

            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v/R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v/R_d - 1) * specific_humidity[i]));
    }
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        ml_parcel_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
}

void Sounding::calc_mlncape() {
    std::vector<double> moist_static_energy;
    std::vector<double> saturated_moist_static_energy;
    for (int i = 0; i < size; i++) {
        moist_static_energy.push_back(c_p * temperature[i] + L_v * specific_humidity[i] + g * height[i]);

        double saturation_mixing_ratio = calc_saturation_mixing_ratio(temperature[i], pressure[i]);
        double qsat = (1 - saturation_mixing_ratio) * saturation_mixing_ratio;
        saturated_moist_static_energy.push_back(c_p * temperature[i] + L_v * qsat + g * height[i]);
    }

    std::vector<double> moist_static_energy_bar;
    moist_static_energy_bar.push_back(moist_static_energy[0]);

    std::vector<double> MSE_trap_sum;
    std::vector<double> height_trap_sum;
    std::vector<double> MSE_times_height;

    for (int i = 0; i < size - 1; i++) {
        MSE_trap_sum.push_back(moist_static_energy[i] + moist_static_energy[i + 1]);
        height_trap_sum.push_back(height[i + 1] - height[i]);
    }

    for (int i = 0; i < MSE_trap_sum.size(); i++) {
        MSE_times_height.push_back(MSE_trap_sum[i] * height_trap_sum[i]);
    }

    for (int i = 1; i < size; i++) {
        double numerator = Helper_Funcs::sumVector(MSE_times_height, 0, i);
        double denominator = 2 * (height[i] - height[0]);
        moist_static_energy_bar.push_back(numerator / denominator);
    }

    
    
    std::vector<double> int_arg;

    for (int i = 0; i < size; i++) {
        double term_one = g / (c_p * temperature[i]);
        double term_two = moist_static_energy_bar[i] - saturated_moist_static_energy[i];
        int_arg.push_back(-1.0 * term_one * term_two);
    }
    

    double sum = 0.0;
    int j = 0;
    for (int i = 0; height[i+1] < ml_el; i++) {
        if(height[i] < ml_lfc) {
            if(height[i+1] < ml_lfc) {
                continue;
            }
            double height_diff = (height[i+1] - ml_lfc);
            double height_prop = height_diff / (height[i+1] - height[i]);
            double int_arg_diff = int_arg[i+1] - int_arg[i];
            double interim_int_arg = (int_arg_diff * height_prop) + int_arg[i];
            double term_one = 0.5 * (interim_int_arg + int_arg[i + 1]);
            double term_two = height[i + 1] - ml_lfc;
            sum += term_one * term_two;
            continue;
            
        }
        double term_one = 0.5 * (int_arg[i] + int_arg[i + 1]);
        double term_two = height[i + 1] - height[i];
        sum += term_one * term_two;
        j = i;
    }
    j++;
    double height_diff = (height[j+1] - ml_el);
    double height_prop = height_diff / (height[j+1] - height[j]);
    double int_arg_diff = int_arg[j+1] - int_arg[j];
    double interim_int_arg = (int_arg_diff * height_prop) + int_arg[j];
    double term_one = 0.5 * (int_arg[j] + interim_int_arg);
    double term_two = ml_el - height[j];
    sum += term_one * term_two;
    
    if(sum < 0) {
        sum = 0;
    }
    
    ml_ncape = sum;
}

void Sounding::calc_ml_e_tilde_and_frac_entrainment() {
    double v_sr = get_srw(0, 1000);
    double H = ml_el;
    double sigma = 1.1;
    double alpha = 0.8;
    double prandtl_number = 1.0 / 3.0;
    double von_karman_constant = 0.18;
    double psi = von_karman_constant * (alpha * alpha) * (std::numbers::pi * std::numbers::pi) * 120.0 / (4.0 * prandtl_number * (sigma * sigma) * H);
    double vsr_tilde = v_sr / (std::pow(2.0 * ml_cape, 0.5));
    double n_tilde = ml_ncape / ml_cape;
    
    double vsr_tilde_squared = vsr_tilde * vsr_tilde;

    double numerator = -1 - psi - (psi / vsr_tilde_squared) * n_tilde;
    double root_term_one = 1 + psi + (psi / vsr_tilde_squared) * n_tilde;
    double root_term_two = 4 * psi / vsr_tilde_squared * (1 - psi * n_tilde);
    double square_root = std::pow(root_term_one * root_term_one + root_term_two, 0.5);

    numerator += square_root;

    double fraction = numerator / (2 * psi / vsr_tilde_squared);
    ml_e_tilde = fraction;

    double e_tilde = fraction;

    ml_fractional_entrainment = 2 * ((1 - e_tilde) / (e_tilde + n_tilde)) / (ml_el);
    
    ml_e_a = vsr_tilde_squared + (-1 - psi - (psi / vsr_tilde_squared) * n_tilde + std::pow(std::pow(1 + psi + (psi / vsr_tilde_squared) * n_tilde, 2) + 4 * (psi / vsr_tilde_squared) * (1 - psi * n_tilde), 0.5)) / (2 * psi / (vsr_tilde_squared));
}

void Sounding::calc_ml_entrainment_buoyancy() {
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    T_lif.push_back(mixed_parcel_start_temp);
    double starting_mixing_ratio = calc_saturation_mixing_ratio(mixed_parcel_start_dew, mixed_parcel_start_pressure);
    double starting_specific_humidity = starting_mixing_ratio / (1 + starting_mixing_ratio);
    Q_v_lif.push_back(starting_specific_humidity);
    Q_t_lif.push_back(starting_specific_humidity);

    double q_sat_prev = 0;
    double i = 0;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if ((Q_v_lif[i - 1] < q_sat) && (pressure[i] > 50000)){
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], ml_fractional_entrainment));
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * ml_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            if (Q_v_lif[i] >= q_sat) {
                
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                
                if(satrat < 0.05) {
                    q_sat_prev = q_sat;
                    continue;
                }
                
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);


                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], ml_fractional_entrainment);
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * ml_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + specific_humidity[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep, ml_fractional_entrainment);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * ml_fractional_entrainment * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1], ml_fractional_entrainment));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (ml_fractional_entrainment * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));
            
            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }
    
    

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v / R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v / R_d - 1) * specific_humidity[i]));
    }
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        ml_entrainment_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
}



void Sounding::calc_muparcel_path() {
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    for(int i = 0; height[i] <= most_unstable_parcel_start_height; i++) {
        T_lif.push_back(temperature[i]);
        Q_v_lif.push_back(specific_humidity[i]);
        Q_t_lif.push_back(specific_humidity[i]);
    }

    double q_sat_prev = 0;
    double i = most_unstable_parcel_start_index;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if (Q_v_lif[i-1]<q_sat) {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0));
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            
            if (Q_v_lif[i] >= q_sat) {
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);



                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], 0);
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * 0 * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + pressure[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters_no_e(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * 0 * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters_no_e(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1]));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (0 * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));

            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v/R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v/R_d - 1) * specific_humidity[i]));
    }
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        mu_parcel_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
}

void Sounding::calc_muncape() {
    std::vector<double> moist_static_energy;
    std::vector<double> saturated_moist_static_energy;
    for (int i = 0; i < size; i++) {
        moist_static_energy.push_back(c_p * temperature[i] + L_v * specific_humidity[i] + g * height[i]);

        double saturation_mixing_ratio = calc_saturation_mixing_ratio(temperature[i], pressure[i]);
        double qsat = (1 - saturation_mixing_ratio) * saturation_mixing_ratio;
        saturated_moist_static_energy.push_back(c_p * temperature[i] + L_v * qsat + g * height[i]);
    }

    std::vector<double> moist_static_energy_bar;
    moist_static_energy_bar.push_back(moist_static_energy[0]);

    std::vector<double> MSE_trap_sum;
    std::vector<double> height_trap_sum;
    std::vector<double> MSE_times_height;

    for (int i = 0; i < size - 1; i++) {
        MSE_trap_sum.push_back(moist_static_energy[i] + moist_static_energy[i + 1]);
        height_trap_sum.push_back(height[i + 1] - height[i]);
    }

    for (int i = 0; i < MSE_trap_sum.size(); i++) {
        MSE_times_height.push_back(MSE_trap_sum[i] * height_trap_sum[i]);
    }

    for (int i = 1; i < size; i++) {
        double numerator = Helper_Funcs::sumVector(MSE_times_height, 0, i);
        double denominator = 2 * (height[i] - height[0]);
        moist_static_energy_bar.push_back(numerator / denominator);
    }

    
    
    std::vector<double> int_arg;

    for (int i = 0; i < size; i++) {
        double term_one = g / (c_p * temperature[i]);
        double term_two = moist_static_energy_bar[i] - saturated_moist_static_energy[i];
        int_arg.push_back(-1.0 * term_one * term_two);
    }
    

    double sum = 0.0;
    int j = 0;
    for (int i = 0; height[i+1] < mu_el; i++) {
        if(height[i] < mu_lfc) {
            if(height[i+1] < mu_lfc) {
                continue;
            }
            double height_diff = (height[i+1] - mu_lfc);
            double height_prop = height_diff / (height[i+1] - height[i]);
            double int_arg_diff = int_arg[i+1] - int_arg[i];
            double interim_int_arg = (int_arg_diff * height_prop) + int_arg[i];
            double term_one = 0.5 * (interim_int_arg + int_arg[i + 1]);
            double term_two = height[i + 1] - mu_lfc;
            sum += term_one * term_two;
            continue;
            
        }
        double term_one = 0.5 * (int_arg[i] + int_arg[i + 1]);
        double term_two = height[i + 1] - height[i];
        sum += term_one * term_two;
        j = i;
    }
    j++;
    double height_diff = (height[j+1] - mu_el);
    double height_prop = height_diff / (height[j+1] - height[j]);
    double int_arg_diff = int_arg[j+1] - int_arg[j];
    double interim_int_arg = (int_arg_diff * height_prop) + int_arg[j];
    double term_one = 0.5 * (int_arg[j] + interim_int_arg);
    double term_two = mu_el - height[j];
    sum += term_one * term_two;
    
    if(sum < 0) {
        sum = 0;
    }
    
    mu_ncape = sum;
}

void Sounding::calc_mu_e_tilde_and_frac_entrainment() {
    double v_sr = get_srw(0, 1000);
    double H = mu_el;
    double sigma = 1.1;
    double alpha = 0.8;
    double prandtl_number = 1.0 / 3.0;
    double von_karman_constant = 0.18;
    double psi = von_karman_constant * (alpha * alpha) * (std::numbers::pi * std::numbers::pi) * 120.0 / (4.0 * prandtl_number * (sigma * sigma) * H);
    double vsr_tilde = v_sr / (std::pow(2.0 * mu_cape, 0.5));
    double n_tilde = mu_ncape / mu_cape;
    
    double vsr_tilde_squared = vsr_tilde * vsr_tilde;

    double numerator = -1 - psi - (psi / vsr_tilde_squared) * n_tilde;
    double root_term_one = 1 + psi + (psi / vsr_tilde_squared) * n_tilde;
    double root_term_two = 4 * psi / vsr_tilde_squared * (1 - psi * n_tilde);
    double square_root = std::pow(root_term_one * root_term_one + root_term_two, 0.5);

    numerator += square_root;

    double fraction = numerator / (2 * psi / vsr_tilde_squared);
    mu_e_tilde =  fraction;

    double e_tilde = fraction;

    mu_fractional_entrainment = 2 * ((1 - e_tilde) / (e_tilde + n_tilde)) / (mu_el);
    
    mu_e_a = vsr_tilde_squared + (-1 - psi - (psi / vsr_tilde_squared) * n_tilde + std::pow(std::pow(1 + psi + (psi / vsr_tilde_squared) * n_tilde, 2) + 4 * (psi / vsr_tilde_squared) * (1 - psi * n_tilde), 0.5)) / (2 * psi / (vsr_tilde_squared));
}

void Sounding::calc_mu_entrainment_buoyancy() {
    std::vector<double> T_lif;
    std::vector<double> Q_v_lif;
    std::vector<double> Q_t_lif;
    for(int i = 0; height[i] <= most_unstable_parcel_start_height; i++) {
        T_lif.push_back(temperature[i]);
        Q_v_lif.push_back(specific_humidity[i]);
        Q_t_lif.push_back(specific_humidity[i]);
    }

    double q_sat_prev = 0;
    double i = most_unstable_parcel_start_index;

    while (i < size - 1) {
        i++;
        double q_sat = (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i - 1], pressure[i - 1]);

        if ((Q_v_lif[i - 1] < q_sat) && (pressure[i] > 50000)){
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], mu_fractional_entrainment));
            Q_v_lif.push_back(Q_v_lif[i - 1] - (height[i] - height[i - 1]) * mu_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]));
            Q_t_lif.push_back(Q_v_lif[i]);
            q_sat = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);
            
            if (Q_v_lif[i] >= q_sat) {
                
                double satrat = (Q_v_lif[i] - q_sat_prev) / (q_sat - q_sat_prev);
                
                if(satrat < 0.05) {
                    q_sat_prev = q_sat;
                    continue;
                }
                
                double dz_dry = satrat * (height[i] - height[i - 1]);
                double dz_wet = (1 - satrat) * (height[i] - height[i - 1]);


                double T_halfstep = T_lif[i - 1] + dz_dry * dry_lapse_rate(T_lif[i - 1], Q_v_lif[i - 1], temperature[i - 1], specific_humidity[i - 1], mu_fractional_entrainment);
                double Qv_halfstep = Q_v_lif[i - 1] - dz_dry * mu_fractional_entrainment * (Q_v_lif[i - 1] - specific_humidity[i - 1]);
                double Qt_halfstep = Q_v_lif[i];
                double p_halfstep = pressure[i - 1] * satrat + specific_humidity[i] * (1 - satrat);
                double T0_halfstep = temperature[i - 1] * satrat + temperature[i] * (1 - satrat);
                double Q0_halfstep = specific_humidity[i - 1] * satrat + specific_humidity[i] * (1 - satrat);

                T_lif[i] = T_halfstep + dz_wet * calc_moist_lapse_rate_peters(T_halfstep, Qv_halfstep, (1 - Qt_halfstep) * compute_mixing_ratio_saturation_linear_combo(T_halfstep, p_halfstep), (1 - Qt_halfstep) * compute_mixing_ratio_saturation_all_ice(T_halfstep, p_halfstep), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep, mu_fractional_entrainment);


                Q_t_lif[i] = Q_t_lif[i - 1] - (height[i] - height[i - 1]) * mu_fractional_entrainment * (Qt_halfstep - Q0_halfstep);
                Q_v_lif[i] = (1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]);

                if (Q_t_lif[i] < Q_v_lif[i]) {
                    Q_v_lif[i] = Q_t_lif[i];
                }
            }
            q_sat_prev = q_sat;
        }
        else {
            T_lif.push_back(T_lif[i - 1] + (height[i] - height[i - 1]) * calc_moist_lapse_rate_peters(T_lif[i - 1], Q_v_lif[i - 1], (1 - Q_t_lif[i - 1]) * calc_saturation_mixing_ratio(T_lif[i - 1], pressure[i - 1]), (1 - Q_t_lif[i - 1]) * compute_mixing_ratio_saturation_all_ice(T_lif[i - 1], pressure[i - 1]), pressure[i - 1], temperature[i - 1], specific_humidity[i - 1], Q_t_lif[i - 1], mu_fractional_entrainment));
            Q_t_lif.push_back(Q_t_lif[i - 1] - (height[i] - height[i - 1]) * (mu_fractional_entrainment * (Q_t_lif[i - 1] - specific_humidity[i - 1])));
            Q_v_lif.push_back((1 - Q_t_lif[i]) * compute_mixing_ratio_saturation_linear_combo(T_lif[i], pressure[i]));
            
            if (Q_t_lif[i] < Q_v_lif[i]) {
                Q_v_lif[i] = Q_t_lif[i];
            }
        }
    }
    
    

    std::vector<double> T_rho_lif;
    std::vector<double> T_0_lif;
    for (int i = 0; i < T_lif.size(); i++) {
        T_rho_lif.push_back(T_lif[i] * (1.0 + (R_v / R_d) * Q_v_lif[i] - Q_t_lif[i]));
        T_0_lif.push_back(temperature[i] * (1.0 + (R_v / R_d - 1) * specific_humidity[i]));
    }
    
    for (int i = 0; i < T_rho_lif.size(); i++) {
        mu_entrainment_buoyancy.push_back(g * (T_rho_lif[i] - T_0_lif[i]) / T_0_lif[i]);
    }
}
