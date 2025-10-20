/**
 * @file butterworth_design.hpp
 * @brief Defines Butterworth Filter implementation with low-pass, high-pass and band-pass types of operation
 *
 * @see https://www.mathworks.com/help/signal/ref/butter.html
 * 
 * @author Gabriel Cardoso da Silva
 * @date Octuber 7, 2025
 */

#ifndef BUTTERWORTH_DESIGN_HPP
#define BUTTERWORTH_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
	template<typename precision, size_t order>
	class ButterworthLowpass : public IIRGenericDesign<precision, order, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* Butterworth low-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Cutoff frequency for the low-pass filter */
		ButterworthLowpass(const precision& sample_rate, const precision& first_cut) {
			//Assigns the sampling frequency
			this->m_sample_rate = sample_rate;

			//Pre-warps the frequencies to the analog domain
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);
			
			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calculates all poles and their conjugates for the low-pass filter prototype
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Stores the calculated pole
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Stores the conjugate right after it
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//If there is a real pole (in case of odd order)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Determines the analog cutoff frequency
			const precision omega = this->m_first_cutoff;

			//The number of analog poles equals the filter order
			const size_t num_analogue_poles = order;

			//Calculates the digital poles by multiplying them by the cutoff frequency
			for (size_t index = 0; index < num_analogue_poles; index++) {
				this->m_poles[index] *= omega;

				//Checks if all poles are in the left half-plane
				if (this->m_poles[index].real() > 0) {
					//Otherwise, stops the process and throws an exception
					return false;
				}
			}

			//Calculates the initial gain of the filter
			precision initial_gain = pow(omega, num_analogue_poles);

			//Performs the plane conversion, calculates the zeros, and obtains the gain after conversion
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Applies a gain correction
			this->m_overall_gain = initial_gain * (initial_gain / this->m_overall_gain);

			return true;
		}

		auto createSOS() -> bool override {
			//Generates the second-order sections from the poles and zeros
			zpk2Biquads(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_biquads, this->m_num_biquads
			);

			//If it reached this point, indicates success
			return true;
		}

	};

	template<typename precision, size_t order>
	class ButterworthHighpass : public IIRGenericDesign<precision, order, order, (order / 2) + (order % 2 == 1)> {
	public:
		/* Butterworth high-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Cutoff frequency for the high-pass filter */
		ButterworthHighpass(const precision& sample_rate, const precision& first_cut) {
			//Assigns the sampling frequency
			this->m_sample_rate = sample_rate;

			//Pre-warps the frequencies to the analog domain
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calculates all poles and their conjugates for the low-pass filter prototype
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Stores the calculated pole
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Stores the conjugate right after it
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//If there is a real pole (in case of odd order)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Determines the analog cutoff frequency
			const precision omega = this->m_first_cutoff;

			//The number of analog poles equals the filter order
			const size_t num_analogue_poles = order;

			//Calculates the digital poles by multiplying them by the cutoff frequency
			for (size_t index = 0; index < num_analogue_poles; index++) {
				//Ignore null poles
				if (cabs(this->m_poles[index]) == 0) continue;

				//Converts the low-pass pole to high-pass using the cutoff frequency
				this->m_poles[index] = Complex<precision>(omega) / this->m_poles[index];

				//Initializes the zeros
				this->m_zeros[index] = Complex<precision>(0.0);

				//Checks if all poles are in the left half-plane
				if (this->m_poles[index].real() > 0) {
					//Otherwise, stops the process and throws an exception
					return false;
				}
			}

			//Calculates the initial gain of the filter
			precision initial_gain = 1.0;

			//Performs the plane conversion, calculates the zeros, and obtains the gain after conversion
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Applies a gain correction
			this->m_overall_gain = precision(1) / this->m_overall_gain;

			return true;
		}

		auto createSOS() -> bool override {
			//Generates the second-order sections from the poles and zeros
			zpk2Biquads<precision>(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_biquads, this->m_num_biquads
			);

			//If it reached this point, indicates success
			return true;
		}

	};

	template<typename precision, size_t order>
	class ButterworthBandpass : public IIRGenericDesign<precision, order * 2, order, order> {
	public:
		/* Butterworth band-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Lower cutoff frequency for the band-pass filter
		 * @param precision& [Hz] Upper cutoff frequency for the band-pass filter */
		ButterworthBandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut) {
			//Assigns the sampling frequency
			this->m_sample_rate = sample_rate;
			
			//Pre-warps the frequencies to the analog domain
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);
			this->m_second_cutoff = 2 * tan(M_PI * second_cut / sample_rate);
			
			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calculates all poles and their conjugates for the low-pass filter prototype
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Stores the calculated pole
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Stores the conjugate right after it
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//If there is a real pole (in case of odd order)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Calculates the bandwidth and the analog cutoff frequency
			const precision bandwidth = this->m_second_cutoff - this->m_first_cutoff;
			const precision omega = sqrt(this->m_first_cutoff * this->m_second_cutoff);

			//The number of analog poles equals the filter order
			const size_t num_analogue_poles = order;

			//Transforms the low-pass prototype into a band-pass filter
			//The number of poles in the band-pass filter is 2n
			for (size_t index = 0; index < num_analogue_poles; index++) {
				//Ignores null poles
				if (cabs(this->m_poles[index]) == 0) continue;

				//If it is a real pole (last pole in the odd-order filter)
				if (this->m_poles[index].imag() == 0 && index == num_analogue_poles - 1 && (order % 2 == 1)) {
					Complex<precision> first_half  = precision(0.5) * this->m_poles[index] * bandwidth + Complex<precision>(0, omega);
					Complex<precision> second_half = precision(0.5) * this->m_poles[index] * bandwidth - Complex<precision>(0, omega);

					this->m_poles[index] = first_half;
					this->m_poles[index + num_analogue_poles] = second_half;

					//Checks if all poles are in the left half-plane
					if (first_half.real() > 0 || second_half.real() > 0) {
						//Otherwise, stops the process and throws an exception
						return false;
					}
				}

				//If it is a complex pole
				else {
					//Calculates the first and second terms
					Complex<precision> first  = precision(0.5) * this->m_poles[index] * bandwidth;
					Complex<precision> second = precision(0.5) * csqrt((bandwidth * bandwidth) * (this->m_poles[index] * this->m_poles[index]) - 4 * omega * omega);

					//The first half should be the product of the first and second terms
					const auto& first_half_pole = this->m_poles[index] = first + second;

					//The second half should be the difference between the first and second terms
					const auto& second_half_pole = this->m_poles[index + num_analogue_poles] = first - second;

					//Checks if all poles are in the left half-plane
					if (first_half_pole.real() > 0 || second_half_pole.real() > 0) {
						//Otherwise, stops the process and throws an exception
						return false;
					}
				}
			}

			//Calculates the initial gain of the filter
			precision initial_gain = pow(bandwidth, num_analogue_poles);

			//Performs the plane conversion, calculates the zeros, and obtains the gain after conversion
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Applies a gain correction
			this->m_overall_gain = initial_gain * (initial_gain / this->m_overall_gain);

			return true;
		}

		auto createSOS() -> bool override {
			//Generates the second-order sections from the poles and zeros
			zpk2Biquads<precision>(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_biquads, this->m_num_biquads
			);

			//If it reached this point, indicates success
			return true;
		}

	};
}
#endif
