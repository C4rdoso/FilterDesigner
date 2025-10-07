/**
 * @file math_utils.hpp
 * @brief Define some math tools tha are required in this project
 *
 * @author Gabriel Cardoso da Silva
 * @date Octuber 7, 2025
 */

#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

//Include the platform definitions
#include "platform_definitions.hpp"

//If the system has support with STL
#if HAS_STL_SUPPORT
#	include <cmath>
#	include <cstdint>
#	include <complex>
#endif

namespace FilterDesigner {

#if !HAS_STL_SUPPORT
	template <typename precision>
	struct Complex {
		//Real and imaginary parts
		precision m_real;
		precision m_imag;

		/* Default constructor for complex number
		 * @param precision Real part of the complex number 
		 * @param precision Imaginary part of the complex number */
		Complex(const precision& real = 0, const precision& imag = 0) : m_real(real), m_imag(imag) {}


		auto operator+(const Complex& other) const -> Complex {
			return Complex(m_real + other.m_real, m_imag + other.m_imag);
		}
		
		auto operator-(const Complex& other) const -> Complex {
			return Complex(m_real - other.m_real, m_imag - other.m_imag);
		}

		auto operator*(const Complex& other) const -> Complex {
			return Complex(
				m_real * other.m_real - m_imag * other.m_imag,
				m_real * other.m_imag + m_imag * other.m_real
			);
		}

		auto operator/(const Complex& other) const -> Complex {
			precision denom = other.m_real * other.m_real + other.m_imag * other.m_imag;
			return Complex(
				(m_real * other.m_real + m_imag * other.m_imag) / denom,
				(m_imag * other.m_real - m_real * other.m_imag) / denom
			);
		}

		// --- Complex and real number operations --- 
		auto operator+(precision value) const -> Complex {
			return Complex(m_real + value, m_imag);
		}
		auto operator-(precision value) const -> Complex {
			return Complex(m_real - value, m_imag);
		}
		auto operator*(precision value) const -> Complex {
			return Complex(m_real * value, m_imag * value);
		}
		auto operator/(precision value) const -> Complex {
			return Complex(m_real / value, m_imag / value);
		}

		// --- Compound assignments (+=, -=, *=, /=) --- 
		auto operator+=(const Complex& other) -> Complex& {
			m_real += other.m_real; m_imag += other.m_imag; return *this;
		}
		auto operator-=(const Complex& other) -> Complex& {
			m_real -= other.m_real; m_imag -= other.m_imag; return *this;
		}
		auto operator*=(const Complex& other) -> Complex& {
			precision r = m_real * other.m_real - m_imag * other.m_imag;
			precision i = m_real * other.m_imag + m_imag * other.m_real;
			m_real = r; m_imag = i;
			return *this;
		}
		auto operator/=(const Complex& other) -> Complex& {
			precision denom = other.m_real * other.m_real + other.m_imag * other.m_imag;
			precision r = (m_real * other.m_real + m_imag * other.m_imag) / denom;
			precision i = (m_imag * other.m_real - m_real * other.m_imag) / denom;
			m_real = r; m_imag = i;
			return *this;
		}

		// --- Compound assignments between complex and real number --- 
		auto operator+=(precision value) -> Complex& { m_real += value; return *this; }
		auto operator-=(precision value) -> Complex& { m_real -= value; return *this; }
		auto operator*=(precision value) -> Complex& { m_real *= value; m_imag *= value; return *this; }
		auto operator/=(precision value) -> Complex& { m_real /= value; m_imag /= value; return *this; }
		
		/* @return The value of the real part */
		auto real() const -> precision { return m_real; }
		/* @return The value of the imaginary part */
		auto imag() const -> precision { return m_imag; }
	};


	// --- Global operations between real and complex number (real number on the left side) ---
	template <typename precision>
	auto operator+(precision lhs, const Complex<precision>& rhs) -> Complex<precision> {
		return Complex<precision>(lhs + rhs.m_real, rhs.m_imag);
	}

	template <typename precision>
	auto operator-(precision lhs, const Complex<precision>& rhs) -> Complex<precision> {
		return Complex<precision>(lhs - rhs.m_real, -rhs.m_imag);
	}

	template <typename precision>
	auto operator*(precision lhs, const Complex<precision>& rhs) -> Complex<precision> {
		return Complex<precision>(lhs * rhs.m_real, lhs * rhs.m_imag);
	}

	template <typename precision>
	auto operator/(precision lhs, const Complex<precision>& rhs) -> Complex<precision> {
		precision denom = rhs.m_real * rhs.m_real + rhs.m_imag * rhs.m_imag;
		return Complex<precision>((lhs * rhs.m_real) / denom, (-lhs * rhs.m_imag) / denom);
	}

	/* Magnetude of a complex number
	 * @param Complex<precision>& Complex number
	 * @return The magnetude of the complex number */
	template<typename precision>
	auto cabs(const Complex<precision>& target) -> precision { 
		return sqrt(target.m_real * target.m_real + target.m_imag * target.m_imag); 
	}

	/* Square root for complex numbers
	 * @param Complex<precision>& Complex number
	 * @return Square root result */
	template<typename precision>
	auto csqrt(const Complex<precision>& value) -> Complex<precision> {
		precision magnitude = sqrt(sqrt(value.m_real * value.m_real + value.m_imag * value.m_imag));
		precision angle = atan2(value.m_imag, value.m_real) / precision(2);
		return Complex<precision>(magnitude * cos(angle), magnitude * sin(angle));
	}

#else
	//Use STL version of the complex number
	template<typename precision>
	using Complex = std::complex<precision>;

	//Use STL max function for complex and real numbers
	using std::max;

	/* Magnetude of a complex number
	 * @param Complex<precision>& Complex number
	 * @return The magnetude of the complex number */
	template<typename precision>
	auto cabs(const Complex<precision>& target) -> precision { 
		return std::abs(target);
	}

	/* Square root for complex numbers
	 * @param Complex<precision>& Complex number
	 * @return Square root result */
	template<typename precision>
	auto csqrt(const Complex<precision>& value) -> Complex<precision> {
		precision magnitude = std::sqrt(std::sqrt(value.real() * value.real() + value.imag() * value.imag()));
		precision angle = atan2(value.imag(), value.real()) / precision(2);
		return Complex<precision>(magnitude * cos(angle), magnitude * sin(angle));
	}

#endif

	/* DF2 Biquad class defintion for real time signal processing
	 * @see https://www.mathworks.com/help/dsphdl/ref/biquadfilter.html */
	template<typename precision>
	struct BiquadSection {
		//Coefficients
		precision m_b0, m_b1, m_b2, m_a1, m_a2;

		//Time memory
		precision m_z1 = 0.0, m_z2 = 0.0;

		/* Default constructor for DF2 Biquad */
		BiquadSection()
		: m_b0(0.0), m_b1(0.0), m_b2(0.0), m_a1(0.0), m_a2(0.0) { }

		/* Creates a DF2 Biquad
		 * @param precision& b Numerator coefficients
		 * @param precision& a Denominator coefficients */
		BiquadSection(const precision& b0, const precision& b1, const precision& b2, const precision& a0, const precision& a1, const precision& a2)
		: m_b0(b0 / a0), m_b1(b1 / a0), m_b2(b2 / a0), m_a1((-a1) / a0), m_a2((-a2) / a0) { }

		/* Signal process
		 * @param precision& x Input signal
		 * @return The result signal */
		auto process(const precision& x) -> precision {
			precision output = m_b0 * x + m_z1;
			m_z1 = m_b1 * x - m_a1 * output + m_z2;
			m_z2 = m_b2 * x - m_a2 * output;
			return output;
		}
	};

	template<typename precision>
	inline auto bilinearTransform(Complex<precision>& value) -> precision {
		Complex<precision> two(2, 0);
		Complex<precision> s = value;
		value = (two + s) / (two - s);
		return cabs((two - s));
	}

	template<typename precision>
	inline auto planeConversion(const precision& initial_gain, Complex<precision>* poles, const size_t num_poles, Complex<precision>* zeros, const size_t num_zeros) -> precision {
		//The gain after conversion must be an adjustment of the initial gain
		double blt_gain = initial_gain;

		//Adjusts the gain and applies the bilinear transform of the zeros
		for (uint32_t i = 0; i < num_zeros; i++) {
			blt_gain /= bilinearTransform<precision>(zeros[i]);
		}

		//Adjusts the gain and applies the bilinear transform of the poles
		for (uint32_t i = 0; i < num_poles; i++) {
			blt_gain *= bilinearTransform<precision>(poles[i]);
		}

		//Returns the gain after the transformation
		return blt_gain;
	}

	template<typename precision>
	inline auto zpk2Biquads(Complex<precision>* poles, const size_t num_poles, Complex<precision>* zeros, const size_t num_zeros, precision overall_gain, BiquadSection<precision>* biquads, const size_t num_biquads) -> size_t {
		//Filter order
		const size_t filter_order = max(num_zeros, num_poles);

		//Number of required sections
		const size_t num_sections = (filter_order + 1) / 2;

		//Checks if we have enough space in the array
		if (num_sections > num_biquads) return false;

		//Keep track of the last processed section
		size_t last_section_index = 0;

		//Processes the pairs
		for (size_t i = 0; i + 1 < filter_order; i += 2, last_section_index++) {
			const auto first_zero = i < num_zeros ? zeros[i] :			Complex<precision>(-1, 0);
			const auto second_zero = i + 1 < num_zeros ? zeros[i + 1] : Complex<precision>(-1, 0);
			const auto first_pole = i < num_poles ? poles[i] :			Complex<precision>(0, 0);
			const auto second_pole = i + 1 < num_poles ? poles[i + 1] : Complex<precision>(0, 0);

			//Numerator
			precision b0 = 1.0;
			precision b1 = -(first_zero + second_zero).real();
			precision b2 = (first_zero * second_zero).real();

			//Denominator
			precision a0 = 1.0;
			precision a1 = -(first_pole + second_pole).real();
			precision a2 = (first_pole * second_pole).real();

			//Applies the global gain only in the first section
			if (last_section_index == 0) {
				b0 *= overall_gain;
				b1 *= overall_gain;
				b2 *= overall_gain;
			}

			biquads[last_section_index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);
		}

		//If the filter order is odd, the last pole and zero are handled separately
		if (filter_order % 2 == 1) {
			const auto last_zero = (filter_order - 1 < num_zeros) ? zeros[filter_order - 1] : Complex<precision>(-1, 0);
			const auto last_pole = (filter_order - 1 < num_poles) ? poles[filter_order - 1] : Complex<precision>(0, 0);

			precision b0 = 1.0;
			precision b1 = -last_zero.real();
			precision b2 = 0.0;

			precision a0 = 1.0;
			precision a1 = -last_pole.real();
			precision a2 = 0.0;

			//Applies gain in the first section (if it is exactly this one)
			if (last_section_index == 0) {
				b0 *= overall_gain;
				b1 *= overall_gain;
				b2 *= overall_gain;
			}

			biquads[last_section_index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);
			last_section_index++;
		}

		//Indicates to the user how many sections were actually generated
		return last_section_index; 
	}

}

#endif
