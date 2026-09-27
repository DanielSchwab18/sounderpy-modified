#pragma once

#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <string>
#include "Sounding.hpp"




class EveryLevelBuoyancy {
public:
    EveryLevelBuoyancy();
    ~EveryLevelBuoyancy();
    int get_every_level_buoyancy(){
    /*
    std::string file = "/Users/danielschwab/Desktop/Purdue/Research/Chavas Spring 2025/soundings modified interpolated/2014052206_H_100_38.69_-099.98_KS_KDDC_102_00100010_MS1_0_2293_-0061_1166_049.4_042.7_034.0_0567_0430_18.8_15.1_5.2_7.9_1.36_1.0_2.1_capecalc";

    Sounding s = Sounding(file, false, false);
    */


    std::string path = "/Users/danielschwab/Desktop/Purdue/Research/Chavas Spring 2025/soundings modified interpolated 100m"; // Change to your target folder

    //std::string path = "/Users/danielschwab/Desktop/Purdue/Research/Code Verification/";
    std::ofstream outfile("/Users/danielschwab/Desktop/Purdue/Research/Chavas Summer 2025/Sounding Variable Calculator/Sounding Variable Calculator/buoyancy_at_every_level.csv");
    //std::ofstream outfile("/Users/danielschwab/Desktop/Purdue/Research/Chavas Summer 2025/Sounding Variable Calculator/Sounding Variable Calculator/well_mixed_parcel_test.csv");
    std::ofstream errfile("errors.txt");

    errfile << "ERRORS WILL GO HERE\n";

    std::string cur_file = "";


    outfile << "YEAR,MONTH,DAY,HOUR,EVENT,MAGNITUDE,LAT,LON,STATE,";
    outfile << "SB_CAPE,SB_CIN,SB_NCAPE,SB_ECAPE,SB_ECIN,SB_E_TILDE,SB_E_A,SB_epsilon,SB_LFC,SB_EL,";
    outfile << "ML_CAPE,ML_CIN,ML_NCAPE,ML_ECAPE,ML_ECIN,ML_E_TILDE,ML_E_A,ML_epsilon,ML_LFC_P,ML_E_LFC_P,";
    outfile << "MU_CAPE,MU_CIN,MU_NCAPE,MU_ECAPE,MU_ECIN,MU_E_TILDE,MU_E_A,MU_epsilon,";

    for(int i = 0; i<=10000; i+=100) {
        outfile << "SB_B_" << i << ",";
    }
    for(int i = 0; i<=10000; i+=100) {
        outfile << "ML_B_" << i << ",";
    }
    for(int i = 0; i<=10000; i+=100) {
        outfile << "MU_B_" << i << ",";
    }
    for(int i = 0; i<=10000; i+=100) {
        outfile << "SB_EB_" << i << ",";
    }
    for(int i = 0; i<=10000; i+=100) {
        outfile << "ML_EB_" << i << ",";
    }
    for(int i = 0; i<=10000; i+=100) {
        outfile << "MU_EB_" << i << ",";
    }
    outfile << "\n";


    try {
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            if(entry.path().filename().string() == ".DS_Store") {
                continue;
            }
            if (std::filesystem::is_regular_file(entry.path())) {
                cur_file = entry.path().string();

                std::string filename = entry.path().filename().string();

                std::vector<std::string> result = Helper_Funcs::split(filename, '_');


                Sounding s = Sounding(cur_file, false, false, 100);

                outfile << result[0].substr(2, 2) << "," << result[0].substr(4, 2) << "," << result[0].substr(6, 2) << "," << result[0].substr(8, 2) << ",";
                outfile << result[1] << "," << result[2] << "," << result[3] << "," << result[4] << "," << result[5] << ",";

                outfile << s.get_sbcape() << "," << s.get_sbcin() << "," << s.get_sbncape_peters() << "," << s.get_sbecape() << "," << s.get_sbecin() << "," << s.get_sb_e() << "," << s.get_sb_e_a_tilde() << "," << s.get_sb_fractional_entrainment() << "," << s.get_sb_lfc() << "," << s.get_sb_el() << ",";

                outfile << s.get_mlcape() << "," << s.get_mlcin() << "," << s.get_mlNCAPE() << "," << s.get_mlecape() << "," << s.get_mlecin() << "," << s.get_ml_e() << "," << s.get_ml_e_a_tilde() << "," << s.get_ml_fractional_entrainment() << "," << s.get_ml_lfc_p() << "," << s.get_ml_e_lfc_p() << ",";

                outfile << s.get_mucape() << "," << s.get_mucin() << "," << s.get_muNCAPE() << "," << s.get_muecape() << "," << s.get_muecin() << "," << s.get_mu_e() << "," << s.get_mu_e_a_tilde() << "," << s.get_mu_fractional_entrainment() << ",";

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_sb_B_at_height(i) << ",";
                }

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_ml_B_at_height(i) << ",";
                }

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_mu_B_at_height(i) << ",";
                }

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_sb_EB_at_height(i) << ",";
                }

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_ml_EB_at_height(i) << ",";
                }

                for(int i = 0; i<=10000; i+=100) {
                    outfile << s.get_mu_EB_at_height(i) << ",";
                }

                outfile << "\n";

            }
        }
    }
    catch (const std::exception& e) {
        errfile << "Error: " << e.what() << " From file: " << cur_file << "\n";
    }

    outfile.close();
    errfile.close();

    return 0;
}
};

inline EveryLevelBuoyancy::EveryLevelBuoyancy() = default;

inline EveryLevelBuoyancy::~EveryLevelBuoyancy() = default;
