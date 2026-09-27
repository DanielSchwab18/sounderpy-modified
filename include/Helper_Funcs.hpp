#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>




class Helper_Funcs {
public:
    Helper_Funcs();
    ~Helper_Funcs();

    template <typename T>
    static void printVector(std::vector<T> vec) {
        std::cout << "[" << vec[0] << ", ";
        for (int i = 1; i < vec.size() - 1; i++) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec[vec.size() - 1] << "]" << std::endl;
    }

    static double sumVector(std::vector<double> vec, int start, int end) {
        double sum = 0.0;
        for (int i = start; i < end; i++) {
            sum += vec[i];
        }
        return sum;
    }

    static std::vector<std::string> split(std::string str, char delimiter) {
        std::vector<std::string> tokens;
        std::stringstream ss(str);
        std::string item;

        while (std::getline(ss, item, delimiter)) {
            tokens.push_back(item);
        }

        return tokens;
    }
private:

};

Helper_Funcs::Helper_Funcs() {}

Helper_Funcs::~Helper_Funcs() {}
