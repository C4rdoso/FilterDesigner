#include "butterworth_design.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>

template <typename precision>
void println_biquads(const FilterDesigner::BiquadSection<precision>* pdata, const size_t lenght, std::ostream& os) {
    for (size_t index = 0; index < lenght; index++) {
        const auto& sos = pdata[index];    
        os << index << ": " << "[ " << std::fixed << std::setprecision(6) << sos.m_b0 << ", " << sos.m_b1 << ", " << sos.m_b2 << ", " <<  sos.m_a1 << ", " << sos.m_a2 << " ]\n"; 
    }
}

template <typename precision>
void println_array(const precision* pdata, const size_t lenght, std::ostream& os) {
    os << "[ ";
    for (size_t index = 0; index < lenght; index++)
        os << pdata[index] << " ";
    os << "]\n";
}

template<typename precision, size_t num_poles, size_t num_zeros, size_t num_sos>
auto print_design(FilterDesigner::IIRGenericDesign<precision, num_poles, num_zeros, num_sos>& design, std::ostream& os) {
    os << "Dimensions (Zeros, Poles, SOS): { " << design.m_num_zeros << ", " << design.m_num_poles << ", " << design.m_num_sos << " }\n";
    os << "Poles: "; println_array(design.m_poles, design.m_num_poles, os);
    os << "Zeros: "; println_array(design.m_zeros, design.m_num_zeros, os);
    os << "Biquads:\n"; println_biquads(design.m_sos_sections, design.m_num_sos, os);
}

int main(void) {
    std::chrono::steady_clock::time_point start, stop;

    start = std::chrono::steady_clock::now();
    FilterDesigner::ButterworthBandpass<double, 4> filter(100.0, 0.5, 15.0);
    stop = std::chrono::steady_clock::now();

    print_design(filter, std::cout);
    std::cout << "\n\nProcessing time: " << stop - start << std::endl;

    return EXIT_SUCCESS;
}