/**
 * @file chebyshev2_design.hpp
 * @brief Defines Chebyshev type II Filter implementation with low-pass, high-pass and band-pass types of operation
 *
 * @see https://www.mathworks.com/help/signal/ref/cheby2.html
 * 
 * @author Gabriel Cardoso da Silva
 * @date Octuber 23, 2025
 */

#ifndef CHEBYSHEV2_DESIGN_HPP
#define CHEBYSHEV2_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
    template<typename precision, size_t order>
    class Chebyshev2Lowpass : public IIRGenericDesign<precision, order, order, (order / 2) + (order % 2 == 1)> {
    public:
        /* Chebyshev Type 2 low-pass filter constructor
         * @param sample_rate [Hz] Sampling rate of the input signal
         * @param first_cut   [Hz] Cutoff frequency for the low-pass filter
         * @param stopband_atten_db [dB] Stopband attenuation in decibels (typically 20 to 80 dB) */
        Chebyshev2Lowpass(const precision& sample_rate, const precision& first_cut, const precision& stopband_atten_db = precision(40)) : m_stopband_atten_db(stopband_atten_db) {
            //Validates stopband attenuation parameter
            if (stopband_atten_db <= 0) m_stopband_atten_db = precision(40);
            
            //Assigns the sampling frequency
            this->m_sample_rate = sample_rate;
            
            //Pre-warps the frequency to the analog domain
            this->m_first_cutoff = precision(2) * tan(M_PI * first_cut / sample_rate);
            
            //Designs the filter
            this->design();
        }

    protected:
        auto createPrototype() -> bool override {
            //Converts stopband attenuation to epsilon (inverse of Type 1)
            const precision epsilon = precision(1) / sqrt(pow(precision(10), m_stopband_atten_db / precision(10)) - precision(1));
            
            //Calculates the parameter for pole positioning
            const precision sinh_term = asinh(precision(1) / epsilon) / precision(order);
            const precision sinh_val = sinh(sinh_term);
            const precision cosh_val = cosh(sinh_term);
            
            //Generates Chebyshev Type 1 prototype first (to be inverted)
            Complex<precision> type1_poles[order > 0 ? order : 1];
            
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision real = -sinh_val * sin(theta);
                const precision imag = cosh_val * cos(theta);
                
                type1_poles[k * 2] = Complex<precision>(real, imag);
                type1_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
            }
            
            if (order % 2 == 1) type1_poles[order - 1] = Complex<precision>(-sinh_val, precision(0));
            
            //Type 2: Invert Type 1 poles to get Type 2 poles
            for (size_t index = 0; index < order; index++) {
                const precision pole_mag = cabs(type1_poles[index]);
                if (pole_mag < precision(1e-12)) return false;
                
                this->m_poles[index] = precision(1) / type1_poles[index];
                
                //Verify poles are in left half-plane
                if (this->m_poles[index].real() >= 0) return false;
            }
            
            //Type 2: Generate zeros on the imaginary axis
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision zero_imag = precision(1) / cos(theta);
                
                //Zeros are purely imaginary (on jω axis)
                this->m_zeros[k * 2] = Complex<precision>(0, zero_imag);
                this->m_zeros[(k * 2) + 1] = Complex<precision>(0, -zero_imag);
            }
            
            //For odd order, last zero is at infinity (represented as very large value)
            if (order % 2 == 1) this->m_zeros[order - 1] = Complex<precision>(0, precision(1e10));
            
            return true;
        }

        auto convertPrototype() -> bool override {
            //Determines the analog cutoff frequency
            const precision omega = this->m_first_cutoff;
            
            //Scale both poles and zeros by cutoff frequency
            for (size_t index = 0; index < order; index++) {
                this->m_poles[index] *= omega;
                this->m_zeros[index] *= omega;
                
                //Verify stability
                if (this->m_poles[index].real() >= 0) return false;
            }
            
            //Calculate initial gain for Type 2 (DC gain should be 1.0)
            precision initial_gain = precision(1);
            
            //Performs the plane conversion
            this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);
            
            //Normalize to ensure unity gain at DC (passband)
            this->m_overall_gain = precision(1) / this->m_overall_gain;
            
            return true;
        }

        auto createSOS() -> bool override {
            zpk2Biquads<precision>(
                this->m_poles, this->m_num_poles,
                this->m_zeros, this->m_num_zeros,
                this->m_overall_gain,
                this->m_biquads, this->m_num_biquads
            );
            
            return true;
        }

    private:
        precision m_stopband_atten_db;
    };

    template<typename precision, size_t order>
    class Chebyshev2Highpass : public IIRGenericDesign<precision, order, order, (order / 2) + (order % 2 == 1)> {
    public:
        /* Chebyshev Type 2 high-pass filter constructor
         * @param sample_rate [Hz] Sampling rate of the input signal
         * @param first_cut [Hz] Cutoff frequency for the high-pass filter
         * @param stopband_atten_db [dB] Stopband attenuation in decibels (typically 20 to 80 dB) */
        Chebyshev2Highpass(const precision& sample_rate, const precision& first_cut, const precision& stopband_atten_db = precision(40)) : m_stopband_atten_db(stopband_atten_db) {
            //Validates stopband attenuation parameter
            if (stopband_atten_db <= 0) m_stopband_atten_db = precision(40);
            
            //Assigns the sampling frequency
            this->m_sample_rate = sample_rate;
            
            //Pre-warps the frequency to the analog domain
            this->m_first_cutoff = precision(2) * tan(M_PI * first_cut / sample_rate);
            
            //Designs the filter
            this->design();
        }

    protected:
        auto createPrototype() -> bool override {
            //Converts stopband attenuation to epsilon
            const precision epsilon = precision(1) / sqrt(pow(precision(10), m_stopband_atten_db / precision(10)) - precision(1));
            
            //Calculates the parameter for pole positioning
            const precision sinh_term = asinh(precision(1) / epsilon) / precision(order);
            const precision sinh_val = sinh(sinh_term);
            const precision cosh_val = cosh(sinh_term);
            
            //Generates Chebyshev Type 1 prototype first
            Complex<precision> type1_poles[order > 0 ? order : 1];
            
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision real = -sinh_val * sin(theta);
                const precision imag = cosh_val * cos(theta);
                
                type1_poles[k * 2] = Complex<precision>(real, imag);
                type1_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
            }
            
            if (order % 2 == 1) type1_poles[order - 1] = Complex<precision>(-sinh_val, precision(0));
            
            //Invert Type 1 poles to get Type 2 poles
            for (size_t index = 0; index < order; index++) {
                const precision pole_mag = cabs(type1_poles[index]);
                if (pole_mag < precision(1e-12)) return false;
                
                this->m_poles[index] = precision(1) / type1_poles[index];
                
                if (this->m_poles[index].real() >= 0) return false;
            }
            
            //Generate zeros on the imaginary axis
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision zero_imag = precision(1) / cos(theta);
                
                this->m_zeros[k * 2] = Complex<precision>(0, zero_imag);
                this->m_zeros[(k * 2) + 1] = Complex<precision>(0, -zero_imag);
            }
            
            if (order % 2 == 1) this->m_zeros[order - 1] = Complex<precision>(0, precision(1e10));
            
            return true;
        }

        auto convertPrototype() -> bool override {
            //Determines the analog cutoff frequency
            const precision omega = this->m_first_cutoff;
            
            //Low-pass to high-pass transformation for both poles and zeros
            for (size_t index = 0; index < order; index++) {
                const precision pole_mag = cabs(this->m_poles[index]);
                const precision zero_mag = cabs(this->m_zeros[index]);
                
                if (pole_mag < precision(1e-12)) return false;
                if (zero_mag < precision(1e-12)) continue;
                
                //LP to HP: p_hp = omega / p_lp, z_hp = omega / z_lp
                this->m_poles[index] = Complex<precision>(omega) / this->m_poles[index];
                this->m_zeros[index] = Complex<precision>(omega) / this->m_zeros[index];
                
                if (this->m_poles[index].real() >= 0) return false;
            }
            
            //Calculate initial gain (unity at Nyquist frequency)
            precision initial_gain = precision(1);
            
            //Performs the plane conversion
            this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);
            
            //Normalize gain
            this->m_overall_gain = precision(1) / this->m_overall_gain;
            
            return true;
        }

        auto createSOS() -> bool override {
            zpk2Biquads<precision>(
                this->m_poles, this->m_num_poles,
                this->m_zeros, this->m_num_zeros,
                this->m_overall_gain,
                this->m_biquads, this->m_num_biquads
            );
            
            return true;
        }

    private:
        precision m_stopband_atten_db;
    };
    
    template<typename precision, size_t order>
    class Chebyshev2Bandpass : public IIRGenericDesign<precision, order * 2, order * 2, order> {
    public:
        /* Chebyshev Type 2 band-pass filter constructor
         * @param sample_rate [Hz] Sampling rate of the input signal
         * @param first_cut   [Hz] Lower cutoff frequency for the band-pass filter
         * @param second_cut  [Hz] Upper cutoff frequency for the band-pass filter
         * @param stopband_atten_db [dB] Stopband attenuation in decibels (typically 20 to 80 dB) */
        Chebyshev2Bandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut, const precision& stopband_atten_db = precision(40)) : m_stopband_atten_db(stopband_atten_db) {
            //Validates stopband attenuation parameter
            if (stopband_atten_db <= 0) m_stopband_atten_db = precision(40);
            
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
            //Converts stopband attenuation to epsilon
            const precision epsilon = precision(1) / sqrt(pow(precision(10), m_stopband_atten_db / precision(10)) - precision(1));
            
            //Calculates the parameter for pole positioning
            const precision sinh_term = asinh(precision(1) / epsilon) / precision(order);
            const precision sinh_val = sinh(sinh_term);
            const precision cosh_val = cosh(sinh_term);
            
            //Generates Chebyshev Type 1 prototype first
            Complex<precision> type1_poles[order > 0 ? order : 1];
            
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision real = -sinh_val * sin(theta);
                const precision imag = cosh_val * cos(theta);
                
                type1_poles[k * 2] = Complex<precision>(real, imag);
                type1_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
            }
            
            if (order % 2 == 1) type1_poles[order - 1] = Complex<precision>(-sinh_val, precision(0));
            
            //Invert Type 1 poles to get Type 2 poles
            for (size_t index = 0; index < order; index++) {
                const precision pole_mag = cabs(type1_poles[index]);
                if (pole_mag < precision(1e-12)) return false;
                
                this->m_poles[index] = precision(1) / type1_poles[index];
                
                if (this->m_poles[index].real() >= 0) return false;
            }
            
            //Generate zeros on the imaginary axis
            for (size_t k = 0; k < order / 2; k++) {
                const precision theta = precision(2 * k + 1) * M_PI / precision(2 * order);
                const precision zero_imag = precision(1) / cos(theta);
                
                this->m_zeros[k * 2] = Complex<precision>(0, zero_imag);
                this->m_zeros[(k * 2) + 1] = Complex<precision>(0, -zero_imag);
            }
            
            if (order % 2 == 1) this->m_zeros[order - 1] = Complex<precision>(0, precision(1e10));
            
            return true;
        }

        auto convertPrototype() -> bool override {
            //Calculates the bandwidth and center frequency
            const precision bandwidth = this->m_second_cutoff - this->m_first_cutoff;
            const precision omega = sqrt(this->m_first_cutoff * this->m_second_cutoff);
            
            //The number of analog poles/zeros equals the filter order
            const size_t num_analogue_poles = order;
            
            //Temporary storage for transformed poles and zeros
            Complex<precision> temp_poles[order * 2];
            Complex<precision> temp_zeros[order * 2];
            
            //Transform low-pass prototype to band-pass
            //Each pole/zero becomes two poles/zeros
            for (size_t index = 0; index < num_analogue_poles; index++) {
                const precision pole_magnitude = cabs(this->m_poles[index]);
                const precision zero_magnitude = cabs(this->m_zeros[index]);
                
                if (pole_magnitude < precision(1e-12)) return false;
                
                //Transform pole: LP to BP
                //If it's a real pole
                if (abs(this->m_poles[index].imag()) < precision(1e-12) && index == num_analogue_poles - 1 && (order % 2 == 1)) {
                    
                    const Complex<precision> scaled_pole = this->m_poles[index] * bandwidth;
                    const Complex<precision> omega_j = Complex<precision>(0, omega);
                    
                    temp_poles[index] = precision(0.5) * scaled_pole + omega_j;
                    temp_poles[index + num_analogue_poles] = precision(0.5) * scaled_pole - omega_j;
                    
                    if (temp_poles[index].real() >= 0 || temp_poles[index + num_analogue_poles].real() >= 0) return false;
                }
                //Complex pole
                else {
                    const Complex<precision> scaled_pole = this->m_poles[index] * bandwidth;
                    const Complex<precision> discriminant = (bandwidth * bandwidth) * (this->m_poles[index] * this->m_poles[index]) - precision(4) * omega * omega;
                    
                    const Complex<precision> first  = precision(0.5) * scaled_pole;
                    const Complex<precision> second = precision(0.5) * csqrt(discriminant);
                    
                    temp_poles[index] = first + second;
                    temp_poles[index + num_analogue_poles] = first - second;
                    
                    if (temp_poles[index].real() >= 0 || temp_poles[index + num_analogue_poles].real() >= 0) return false;
                }
                
                //Transform zero: LP to BP
                if (zero_magnitude > precision(1e-12) && zero_magnitude < precision(1e9)) {
                    //If it's a real zero (shouldn't occur for Type 2 imaginary zeros)
                    if (abs(this->m_zeros[index].imag()) < precision(1e-12) && index == num_analogue_poles - 1 && (order % 2 == 1)) {
                        
                        const Complex<precision> scaled_zero = this->m_zeros[index] * bandwidth;
                        const Complex<precision> omega_j = Complex<precision>(0, omega);
                        
                        temp_zeros[index] = precision(0.5) * scaled_zero + omega_j;
                        temp_zeros[index + num_analogue_poles] = precision(0.5) * scaled_zero - omega_j;
                    }

                    //Complex zero (imaginary axis zeros)
                    else {
                        const Complex<precision> scaled_zero = this->m_zeros[index] * bandwidth;
                        const Complex<precision> discriminant = (bandwidth * bandwidth) * (this->m_zeros[index] * this->m_zeros[index]) - precision(4) * omega * omega;
                        
                        const Complex<precision> first  = precision(0.5) * scaled_zero;
                        const Complex<precision> second = precision(0.5) * csqrt(discriminant);
                        
                        temp_zeros[index] = first + second;
                        temp_zeros[index + num_analogue_poles] = first - second;
                    }
                }
                else {
                    //Zero at infinity becomes zeros at +/- j*omega (stopband edges)
                    temp_zeros[index] = Complex<precision>(0, omega);
                    temp_zeros[index + num_analogue_poles] = Complex<precision>(0, -omega);
                }
            }
            
            //Copy transformed poles and zeros back
            for (size_t i = 0; i < order * 2; i++) {
                this->m_poles[i] = temp_poles[i];
                this->m_zeros[i] = temp_zeros[i];
            }
            
            //Calculate initial gain
            precision initial_gain = precision(1);
            
            //Performs the plane conversion
            this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);
            
            //Normalize to ensure proper gain in passband
            this->m_overall_gain = precision(1) / this->m_overall_gain;
            
            return true;
        }

        auto createSOS() -> bool override {
            zpk2Biquads<precision>(
                this->m_poles, this->m_num_poles,
                this->m_zeros, this->m_num_zeros,
                this->m_overall_gain,
                this->m_biquads, this->m_num_biquads
            );
            
            return true;
        }

    private:
        precision m_stopband_atten_db;
    };
}

#endif