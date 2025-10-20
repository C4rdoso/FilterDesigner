/**
 * @file butterworth_design.hpp
 * @brief Defines Butterworth Filter implementation with low-pass, high-pass and band-pass types of operation
 *
 * @see https://www.mathworks.com/help/signal/ref/butter.html
 * 
 * @author Gabriel Cardoso da Silva
 * @date Octuber 20, 2025
 */

#ifndef CHEBYSEV1_DESIGN_HPP
#define CHEBYSEV1_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
    
    template<typename precision, size_t order>
    class Chebyshev1Bandpass : public IIRGenericDesign<precision, order * 2, order, order> {
    public:
        /* Chebyshev Type 1 band-pass filter constructor
         * @param sample_rate [Hz] Sampling rate of the input signal
         * @param first_cut   [Hz] Lower cutoff frequency for the band-pass filter
         * @param second_cut  [Hz] Upper cutoff frequency for the band-pass filter
         * @param ripple_db   [dB] Passband ripple in decibels (typically 0.1 to 3 dB) */
        Chebyshev1Bandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut,const precision& ripple_db = precision(0.5)) : m_ripple_db(ripple_db) {
            //Validates ripple parameter, if not valid, assigns a default value
            if (ripple_db <= 0) m_ripple_db = precision(0.5);
            
            //Assigns the sampling frequency
            this->m_sample_rate = sample_rate;
            
            //Pre-warps the frequencies to the analog domain
            this->m_first_cutoff = precision(2) * tan(M_PI * first_cut / sample_rate);
            this->m_second_cutoff = precision(2) * tan(M_PI * second_cut / sample_rate);
            
            //Designs the filter
            this->design();
        }

    protected:
        auto createPrototype() -> bool override {
            //Converts ripple from dB to linear scale
            const precision epsilon     = sqrt(pow(precision(10), m_ripple_db / precision(10)) - precision(1));
            
            //Calculates the parameter for pole positioning
            const precision sinh_term   = asinh(precision(1) / epsilon) / precision(order);
            const precision sinh_val    = sinh(sinh_term);
            const precision cosh_val    = cosh(sinh_term);
            
            //Generates Chebyshev Type 1 low-pass prototype poles
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                
                //Chebyshev poles lie on an ellipse (not a circle like Butterworth)
                const precision real = -sinh_val * sin(theta);
                const precision imag = cosh_val * cos(theta);
                
                //Stores the calculated pole
                this->m_poles[k * 2] = Complex<precision>(real, imag);
                
                //Stores the conjugate right after it
                this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
            }
            
            //If there is a real pole (in case of odd order)
            if (order % 2 == 1) {
                this->m_poles[order - 1] = Complex<precision>(-sinh_val, precision(0));
            }
            
            return true;
        }

        auto convertPrototype() -> bool override {
            //Calculates the bandwidth and the center frequency
            const precision bandwidth = this->m_second_cutoff - this->m_first_cutoff;
            const precision omega = sqrt(this->m_first_cutoff * this->m_second_cutoff);
            
            //The number of analog poles equals the filter order
            const size_t num_analogue_poles = order;
            
            //Transforms the low-pass prototype into a band-pass filter
            //The number of poles in the band-pass filter is 2n
            for (size_t index = 0; index < num_analogue_poles; index++) {
                const precision pole_magnitude = cabs(this->m_poles[index]);
                
                //Ignores null poles (should not occur in Chebyshev)
                if (pole_magnitude < precision(1e-12)) continue;
                
                //If it is a real pole (last pole in the odd-order filter)
                if (abs(this->m_poles[index].imag()) < precision(1e-12) && 
                    index == num_analogue_poles - 1 && 
                    (order % 2 == 1)) {
                    
                    const Complex<precision> scaled_pole = this->m_poles[index] * bandwidth;
                    const Complex<precision> omega_j     = Complex<precision>(0, omega);
                    
                    Complex<precision> first_half  = precision(0.5) * scaled_pole + omega_j;
                    Complex<precision> second_half = precision(0.5) * scaled_pole - omega_j;
                    
                    this->m_poles[index] = first_half;
                    this->m_poles[index + num_analogue_poles] = second_half;
                    
                    //Checks if all poles are in the left half-plane
                    if (first_half.real() >= 0 || second_half.real() >= 0) {
                        //Otherwise, stops the process and throws an exception
                        return false;
                    }
                }
                //If it is a complex pole
                else {
                    //Low-pass to band-pass transformation formula
                    const Complex<precision> scaled_pole    = this->m_poles[index] * bandwidth;
                    const Complex<precision> discriminant   = (bandwidth * bandwidth) * (this->m_poles[index] * this->m_poles[index]) - precision(4) * omega * omega;
                    
                    const Complex<precision> first  = precision(0.5) * scaled_pole;
                    const Complex<precision> second = precision(0.5) * csqrt(discriminant);
                    
                    //The first half pole
                    const auto& first_half_pole = this->m_poles[index] = first + second;
                    
                    //The second half pole
                    const auto& second_half_pole = this->m_poles[index + num_analogue_poles] = first - second;
                    
                    //Checks if all poles are in the left half-plane
                    if (first_half_pole.real() >= 0 || second_half_pole.real() >= 0) {
                        return false;
                    }
                }
                
                //Initializes zeros at the origin (characteristic of bandpass filters)
                if (index < order) this->m_zeros[index] = Complex<precision>(0, 0);
            }
            
            // Calculates the initial gain of the filter
            // For Chebyshev Type 1, gain correction depends on filter order
            precision initial_gain = pow(bandwidth, num_analogue_poles);
            
            // Additional gain correction for Chebyshev passband ripple
            if (order % 2 == 0) {
                // For even order, gain at DC is affected by ripple
                const precision ripple_factor = precision(1) / sqrt(precision(1) + pow(precision(10), m_ripple_db / precision(10)) - precision(1));
                initial_gain *= ripple_factor;
            }
            
            // Performs the plane conversion, calculates the zeros, and obtains the gain after conversion
            this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);
            
            // Applies a gain correction
            this->m_overall_gain = initial_gain * (initial_gain / this->m_overall_gain);
            
            return true;
        }

        auto createSOS() -> bool override {
            // Generates the second-order sections from the poles and zeros
            zpk2Biquads<precision>(
                this->m_poles, this->m_num_poles,
                this->m_zeros, this->m_num_zeros,
                this->m_overall_gain,
                this->m_sos_sections, this->m_num_sos
            );
            
            return true;
        }

    private:
        //Passband ripple in decibels
        precision m_ripple_db;  
    };
}

#endif