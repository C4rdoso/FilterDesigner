#ifndef FILTER_DESIGN_HPP
#define FILTER_DESIGN_HPP

#include "math_utils.hpp"

namespace FilterDesigner {
  
	template<typename precision, size_t num_poles, size_t num_zeros, size_t num_sos>
	class IIRGenericDesign {
	public:
		//Armazena as dimensões dos arrays
		const size_t m_num_poles		{ num_poles };
		const size_t m_num_zeros		{ num_zeros };
		const size_t m_num_sos			{ num_sos };

		//Zeros e polos do filtro
		Complex<precision> m_poles[num_poles > 0 ? num_poles : 1];
		Complex<precision> m_zeros[num_zeros > 0 ? num_zeros : 1];

		//Seções de segunda ordem do filtro
		BiquadSection<precision> m_sos_sections[num_sos > 0 ? num_sos : 1];

		//Ganho inicial e ganho global do filtro
		precision m_overall_gain		{ precision(0) };

		//Parâmetros de configuração do filtro
		precision m_first_cutoff		{ precision(0) };
		precision m_second_cutoff		{ precision(0) };
		precision m_sample_rate			{ precision(0) };

		auto design() -> bool {
			//A primeira frequência de corte deve ser maior que zero
			if (m_first_cutoff <= 0) return false;

			//A segunda frequência de corte deve ser maior que a primeira (caso exista)
			if (m_second_cutoff > 0 && m_second_cutoff <= m_first_cutoff) return false;

			//---------------------------------------------------------------
			// Pipeline de projeto do filtro 
			// 
			// 1o passo: Gera o protótipo analógico
			// 2o passo: Converte o projeto analógico em digital
			// 3o passo: Gera as seções de segunda ordem
			//---------------------------------------------------------------
			if (!createPrototype()) return false;
			if (!convertPrototype()) return false;
			if (!createSOS()) return false;

			//Indica que o projeto foi realizado com sucesso
			m_design_ok = true;

			//Retorna o status do projeto
			return m_design_ok;
		}

		/* @return Verdadeiro se o design do filtro ocorreu com sucesso. */
		auto designSuccess()	const -> bool { return m_design_ok; }

	private:
		//Indica se o projeto foi realizado com sucesso
		bool m_design_ok{ false };

	protected:
		/* Método virtual para criação do protótipo do filtro
		 * @return Verdadeiro se a criação do protótipo ocorrer com sucesso */
		virtual auto createPrototype() -> bool = 0;
		/* Método virtual para conversão do protótipo analógico para digital
		 * @return Verdadeiro se a conversão do protótipo ocorrer com sucesso */
		virtual auto convertPrototype() -> bool = 0;
		/* Método virtual para criação das seções de segunda ordem
		 * @return Verdadeiro se a criação das sessoes ocorrer com sucesso */
		virtual auto createSOS() -> bool = 0;
	};

}

#endif
