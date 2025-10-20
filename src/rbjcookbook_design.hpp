/**
 * @file rbjcookbook_design.hpp
 * @brief Defines RBJ-CookBook Filter implementation with low-pass, high-pass and band-pass types of operation
 *
 * @see https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
 * 
 * @author Gabriel Cardoso da Silva
 * @date Octuber 7, 2025
 */

#ifndef RBJCOOKBOOK_DESIGN_HPP
#define RBJCOOKBOOK_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
	
	template<typename precision, size_t order>
	class RBJCookbookLowpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* RBJ Audio-EQ-Cookbook low-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Cutoff frequency for the low-pass filter */
		RBJCookbookLowpass(const precision& sample_rate, const precision& first_cut) {
			//Assigns the sampling frequency and cutoff frequency
			this->m_sample_rate   = sample_rate;
			this->m_first_cutoff  = first_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//The RBJ filter does not use analog poles or zeros
			return true;
		}

		auto convertPrototype() -> bool override {
			//The RBJ filter does not use prototype conversion
			return true;
		}

		auto createSOS() -> bool override {
      		//Pre-warped digital cutoff frequency
      		const precision omega  = precision(2) * M_PI * this->m_first_cutoff / this->m_sample_rate;
      		const precision cos_w0 = cos(omega);
      		const precision sin_w0 = sin(omega);

      		//Fixed quality (Butterworth 2nd order: ~0.707)
      		const precision quality = precision(1) / sqrt(precision(2));
      		const precision alpha   = sin_w0 / (precision(2) * quality);
      
      		//RBJ low-pass coefficients
      		precision a0 = precision(1) + alpha;
  
      		precision b0 = (precision(1) - cos_w0) / (precision(2) * a0);
      		precision b1 = (precision(1) - cos_w0) / a0;
      		precision b2 = (precision(1) - cos_w0) / (precision(2) * a0);
      		precision a1 = (-precision(2) * cos_w0) / a0;
      		precision a2 = (precision(1) - alpha) / a0;
  
      		//Forms the SOS sections
      		for (size_t index = 0; index < this->m_num_biquads; index++)
      		    this->m_biquads[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//If it reached this point, indicates success
			return true;
		}
	};

	template<typename precision, size_t order>
	class RBJCookbookHighpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* RBJ Audio-EQ-Cookbook high-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Cutoff frequency for the high-pass filter */
		RBJCookbookHighpass(const precision& sample_rate, const precision& first_cut) {
			//Assigns the sampling frequency and cutoff frequency
			this->m_sample_rate = sample_rate;
			this->m_first_cutoff = first_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//The RBJ filter does not use analog poles or zeros
			return true;
		}

		auto convertPrototype() -> bool override {
			//The RBJ filter does not use prototype conversion
			return true;
		}

		auto createSOS() -> bool override {
      		//Pre-warped digital cutoff frequency
      		const precision omega  = precision(2) * M_PI * this->m_first_cutoff / this->m_sample_rate;
      		const precision cos_w0 = cos(omega);
      		const precision sin_w0 = sin(omega);
  
      		//Fixed quality (Butterworth 2nd order: ~0.707)
      		const precision quality = precision(1) / sqrt(precision(2));
      		const precision alpha = sin_w0 / (precision(2) * quality);
  
      		//RBJ high-pass coefficients
      		precision a0 = precision(1) + alpha;
  
      		precision b0 = (precision(1) + cos_w0) / (precision(2) * a0);
      		precision b1 = -(precision(1) + cos_w0) / a0;
      		precision b2 = (precision(1) + cos_w0) / (precision(2) * a0);
      		precision a1 = (-precision(2) * cos_w0) / a0;
      		precision a2 = (precision(1) - alpha) / a0;
  
      		//Forms the SOS sections
      		for (size_t index = 0; index < this->m_num_biquads; index++)
      		    this->m_biquads[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//If it reached this point, indicates success
			return true;
		}
	};

	template<typename precision, size_t order>
	class RBJCookbookBandpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* RBJ Audio-EQ-Cookbook band-pass filter constructor
		 * @param precision& [Hz] Sampling rate of the input signal
		 * @param precision& [Hz] Lower cutoff frequency for the band-pass filter
		 * @param precision& [Hz] Upper cutoff frequency for the band-pass filter */
		RBJCookbookBandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut) {
			//Assigns the sampling frequency and cutoff frequency
			this->m_sample_rate = sample_rate;
			this->m_first_cutoff  = first_cut;
			this->m_second_cutoff = second_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//The RBJ filter does not use analog poles or zeros
			return true;
		}

		auto convertPrototype() -> bool override {
			//The RBJ filter does not use prototype conversion
			return true;
		}

		auto createSOS() -> bool override {
			//Calculates the bandwidth and the analog cutoff frequency
			const precision bandwidth = precision(2) * M_PI * (this->m_second_cutoff -  this->m_first_cutoff) / this->m_sample_rate;
			const precision omega 		= precision(2) * M_PI * sqrt(this->m_first_cutoff * this->m_second_cutoff) / this->m_sample_rate;

      		//ω for RBJ (digital) obtained from analog pre-warping
      		const precision alpha  = sin(omega) * sinh(log(precision(2))/precision(2) * bandwidth / omega);
			const precision cos_w0 = cos(omega);

      		//Analytical RBJ coefficients for band-pass
			precision a0 = precision(1) + alpha;

      		precision b0 = alpha / a0;
      		precision b1 = 0;
      		precision b2 = -alpha / a0;
      		precision a1 = -precision(2) * cos_w0 / a0;
      		precision a2 = (precision(1) - alpha) / a0;

			//Forms the biquads from the second-order coefficients
			for (size_t index = 0; index < this->m_num_biquads; index++)
				this->m_biquads[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//If it reached this point, indicates success
			return true;
		}
	};

}

#endif
