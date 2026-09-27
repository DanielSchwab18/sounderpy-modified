#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <string>

#include "include/Sounding.hpp"
#include "include/EveryLevelBuoyancy.hpp"

namespace fs = std::filesystem;

fs::path get_long_path(const fs::path& relative_or_absolute_path) {
    // 1. Resolve relative path to absolute path
    fs::path abs_path = fs::absolute(relative_or_absolute_path);

    // 2. Format string to add \\?\ if on Windows and not already present
    std::wstring path_str = abs_path.wstring();
    const std::wstring prefix = L"\\\\?\\";
    if (path_str.rfind(prefix, 0) != 0) { // Checks if path_str starts with \\?\
        path_str = prefix + path_str;
    }

    return fs::path(path_str);
}

int main()
{
    
    std::string file = "C:/Users/danie/Desktop/Purdue/Research/Chavas Summer 2025/Sounding Variable Calculator/Sounding Variable Calculator/example";
    
    Sounding s = Sounding(file, true, true, 100);
    std::cout << s.get_srh(0, 500) << std::endl;
    std::cout << s.compute_stp() << std::endl;
    std::cout << s.get_mu_lcl() << std::endl;

    /*
    //EveryLevelBuoyancy e = EveryLevelBuoyancy();
    //e.get_every_level_buoyancy();
    //return 0;
    
    std::string path = "C:/Users/danie/RESEARCH_INPUT/100m interpolated soundings/"; // Change to your target folder

    std::ofstream outfile("../Sounding Variable Calculator/output/windows_test_1.csv");
    std::ofstream errfile("../Sounding Variable Calculator/output/errors.txt");

    errfile << "ERRORS WILL GO HERE\n";
    
    std::string cur_file = "";


    outfile << "YEAR,MONTH,DAY,HOUR,EVENT,MAGNITUDE,LAT,LON,STATE,";
    outfile << "SB_CAPE,SB_CIN,SB_NCAPE,SB_ECAPE,SB_ECIN,SB_E_TILDE,SB_E_A,SB_epsilon,SB_LFC,SB_EL,";
    outfile << "ML_CAPE,ML_CIN,ML_NCAPE,ML_ECAPE,ML_ECIN,ML_E_TILDE,ML_E_A,ML_epsilon,ML_LFC_P,ML_E_LFC_P,";
    outfile << "MU_CAPE,MU_CIN,MU_NCAPE,MU_ECAPE,MU_ECIN,MU_E_TILDE,MU_E_A,MU_epsilon,";
    outfile << "SB_3CAPE,SB_3ECAPE,ML_3CAPE,ML_3ECAPE,MU_3CAPE,MU_3ECAPE,";
    outfile << "RH_0_500,RH_0_1000,RH_0_3000,RH_1000_3000,";
    outfile << "W_0_500,W_0_1000,W_0_3000,W_1000_3000,";
    outfile << "SRW_0_500,SRW_0_1000,SRW_0_3000,SRW_1000_3000,";
    outfile << "SHEAR_0_500,SHEAR_0_1000,SHEAR_0_3000,SHEAR_1000_3000,SHEAR_0_6000,";
    outfile << "SRH_0_500,SRH_0_1000,SRH_0_3000,SRH_1000_3000,SRH_0_6000,";
    outfile << "STREAM_VORT_0_500,STREAM_VORT_0_1000,STREAM_VORT_0_3000,STREAM_VORT_1000_3000,STREAM_VORT_0_6000,";
    outfile << "SB_B_500,SB_B_600,SB_B_700,SB_B_800,SB_B_900,SB_B_1000,SB_B_1100,SB_B_1200,SB_B_1300,SB_B_1400,SB_B_1500,SB_B_2000,SB_B_2500,SB_B_3000,SB_B_3500,SB_B_4000,SB_B_4500,SB_B_5000,";
    outfile << "ML_B_500,ML_B_600,ML_B_700,ML_B_800,ML_B_900,ML_B_1000,ML_B_1100,ML_B_1200,ML_B_1300,ML_B_1400,ML_B_1500,ML_B_2000,ML_B_2500,ML_B_3000,ML_B_3500,ML_B_4000,ML_B_4500,ML_B_5000,";
    outfile << "MU_B_500,MU_B_600,MU_B_700,MU_B_800,MU_B_900,MU_B_1000,MU_B_1100,MU_B_1200,MU_B_1300,MU_B_1400,MU_B_1500,MU_B_2000,MU_B_2500,MU_B_3000,MU_B_3500,MU_B_4000,MU_B_4500,MU_B_5000,";
    outfile << "SB_EB_500,SB_EB_600,SB_EB_700,SB_EB_800,SB_EB_900,SB_EB_1000,SB_EB_1100,SB_EB_1200,SB_EB_1300,SB_EB_1400,SB_EB_1500,SB_EB_2000,SB_EB_2500,SB_EB_3000,SB_EB_3500,SB_EB_4000,SB_EB_4500,SB_EB_5000,";
    outfile << "ML_EB_500,ML_EB_600,ML_EB_700,ML_EB_800,ML_EB_900,ML_EB_1000,ML_EB_1100,ML_EB_1200,ML_EB_1300,ML_EB_1400,ML_EB_1500,ML_EB_2000,ML_EB_2500,ML_EB_3000,ML_EB_3500,ML_EB_4000,ML_EB_4500,ML_EB_5000,";
    outfile << "MU_EB_500,MU_EB_600,MU_EB_700,MU_EB_800,MU_EB_900,MU_EB_1000,MU_EB_1100,MU_EB_1200,MU_EB_1300,MU_EB_1400,MU_EB_1500,MU_EB_2000,MU_EB_2500,MU_EB_3000,MU_EB_3500,MU_EB_4000,MU_EB_4500,MU_EB_5000,";
    outfile << "SB_B_LCL,ML_B_LCL,MU_B_LCL,SB_EB_LCL,ML_EB_LCL,MU_EB_LCL,";
    outfile << "SB_LCL,ML_LCL,MU_LCL,SB_E_LCL,ML_E_LCL,MU_E_LCL,";
    outfile << "SB_LCIN,SB_LECIN,ML_LCIN,ML_LECIN,MU_LCIN,MU_LECIN,ML_TOP,SIZE,STP\n";


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
                
                outfile << s.get_sb3cape() << "," << s.get_sb3ecape() << "," << s.get_ml3cape() << "," << s.get_ml3ecape() << "," << s.get_mu3cape() << "," << s.get_mu3ecape() << ",";

                outfile << s.get_rh(0,500) << "," << s.get_rh(0,1000) << "," << s.get_rh(0,3000) << "," << s.get_rh(1000,3000) << ",";

                outfile << s.get_w(0,500) << "," << s.get_w(0,1000) << "," << s.get_w(0,3000) << "," << s.get_w(1000,3000) << ",";

                outfile << s.get_srw(0,500) << "," << s.get_srw(0,1000) << "," << s.get_srw(0,3000) << "," << s.get_srw(1000,3000) << ",";

                outfile << s.get_shear(0,500) << "," << s.get_shear(0,1000) << "," << s.get_shear(0,3000) << "," << s.get_shear(1000,3000) << "," << s.get_shear(0,6000) << ",";

                outfile << s.get_srh(0,500) << "," << s.get_srh(0,1000) << "," << s.get_srh(0,3000) << "," << s.get_srh(1000,3000) << "," << s.get_srh(0,6000) << ",";

                outfile << s.get_streamwise_vorticity(0,500) << "," << s.get_streamwise_vorticity(0,1000) << "," << s.get_streamwise_vorticity(0,3000) << "," << s.get_streamwise_vorticity(1000,3000) << "," << s.get_streamwise_vorticity(0,6000) << ",";
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_sb_B_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_sb_B_to_height(i) << ",";
                }
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_ml_B_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_ml_B_to_height(i) << ",";
                }
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_mu_B_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_mu_B_to_height(i) << ",";
                }
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_sb_EB_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_sb_EB_to_height(i) << ",";
                }
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_ml_EB_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_ml_EB_to_height(i) << ",";
                }
                
                for(int i = 500; i<=1500; i+=100) {
                    outfile << s.get_mu_EB_to_height(i) << ",";
                }
                for(int i = 2000; i<=5000; i+=500) {
                    outfile << s.get_mu_EB_to_height(i) << ",";
                }
                
                outfile << s.get_sb_B_to_height(s.get_nearest_sb_lcl_height()) << "," << s.get_ml_B_to_height(s.get_nearest_ml_lcl_height()) << ",";
                outfile << s.get_mu_B_to_height(s.get_nearest_mu_lcl_height()) << "," << s.get_sb_EB_to_height(s.get_nearest_sb_e_lcl_height()) << ",";
                outfile << s.get_ml_EB_to_height(s.get_nearest_ml_e_lcl_height()) << "," << s.get_mu_EB_to_height(s.get_nearest_mu_e_lcl_height()) << ",";

                outfile << s.get_sb_lcl() << "," << s.get_ml_lcl() << "," << s.get_mu_lcl() << ",";
                outfile << s.get_sb_e_lcl() << "," << s.get_ml_e_lcl() << "," << s.get_mu_lcl() << ",";

                outfile << s.get_sb_lcin() << "," << s.get_sb_lecin() << "," << s.get_ml_lcin() << "," << s.get_ml_lecin() << "," << s.get_mu_lcin() << "," << s.get_mu_lcin() << "," << s.get_ml_top() << ",";
                
                outfile << s.get_size() << "," << s.compute_stp() << "\n";

            }
        }
    }
    catch (const std::exception& e) {
        errfile << "Error: " << e.what() << " From file: " << cur_file << "\n";
    }
    
    outfile.close();
    errfile.close();
    */
    return 0;
}

