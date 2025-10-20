/**
 * @file filter_design.hpp
 * @brief Defines the virtual class that is inherited by all implemented IIR digital filters
 *
 * @author Gabriel Cardoso da Silva
 * @date Octuber 7, 2025
 */

#ifndef FILTER_DESIGN_HPP
#define FILTER_DESIGN_HPP

#include "math_utils.hpp"

namespace FilterDesigner {
	/* IIR Designer declaration
	 * The number of bytes allocated will be based on template parameters */
	template<typename precision, size_t num_poles, size_t num_zeros, size_t num_biquads>
	class IIRGenericDesign {
	public:
		//Stores the dimensions of the arrays
		const size_t m_num_poles		{ num_poles };
		const size_t m_num_zeros		{ num_zeros };
		const size_t m_num_biquads		{ num_biquads };

		//Filter zeros and poles
		Complex<precision> m_poles[num_poles > 0 ? num_poles : 1];
		Complex<precision> m_zeros[num_zeros > 0 ? num_zeros : 1];

		//Second-order sections of the filter
		BiquadSection<precision> m_biquads[num_biquads > 0 ? num_biquads : 1];

		//Initial gain and global gain of the filter
		precision m_overall_gain		{ precision(0) };

		//Filter configuration parameters
		precision m_first_cutoff		{ precision(0) };
		precision m_second_cutoff		{ precision(0) };
		precision m_sample_rate			{ precision(0) };

		auto design() -> bool {
			//The first cutoff frequency must be greater than zero
			if (m_first_cutoff <= 0) return false;

			//The second cutoff frequency must be greater than the first (if it exists)
			if (m_second_cutoff > 0 && m_second_cutoff <= m_first_cutoff) return false;

			//---------------------------------------------------------------
			// Filter design pipeline 
			// 
			// Step 1: Generate the analog prototype
			// Step 2: Convert the analog design to digital
			// Step 3: Generate the second-order sections
			//---------------------------------------------------------------
			if (!createPrototype()) return false;
			if (!convertPrototype()) return false;
			if (!createSOS()) return false;

			//Indicates that the design was completed successfully
			m_design_ok = true;

			//Returns the design status
			return m_design_ok;
		}

		/* @return True if the filter design was successful. */
		auto designSuccess()	const -> bool { return m_design_ok; }

	private:
		//Indicates whether the design was completed successfully
		bool m_design_ok{ false };

	protected:
		/* Virtual method for creating the filter prototype
		 * @return True if the prototype creation was successful */
		virtual auto createPrototype() -> bool = 0;
		/* Virtual method for converting the analog prototype to digital
		 * @return True if the prototype conversion was successful */
		virtual auto convertPrototype() -> bool = 0;
		/* Virtual method for creating the second-order sections
 		 * @return True if the section creation was successful */
		virtual auto createSOS() -> bool = 0;
	};

}

#endif
