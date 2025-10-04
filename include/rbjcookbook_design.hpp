#ifndef RBJCOOKBOOK_DESIGN_HPP
#define RBJCOOKBOOK_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
	
	template<typename precision, size_t order>
	class RBJCookbookLowpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* Construtor do filtro RBJ Audio-EQ-Cookbook do tipo passa-baixa
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte para passa-baixa */
		RBJCookbookLowpass(const precision& sample_rate, const precision& first_cut) {
			//Atribui a frequência de amostragem e a frequência de corte
			this->m_sample_rate   = sample_rate;
			this->m_first_cutoff  = first_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//O filtro RBJ não utiliza polos ou zeros analógicos
			return true;
		}

		auto convertPrototype() -> bool override {
			//O filtro RBJ não utiliza conversão de protótipo
			return true;
		}

		auto createSOS() -> bool override {
      		//Frequência de corte digital pré-distorcida
      		const precision omega  = precision(2) * M_PI * this->m_first_cutoff / this->m_sample_rate;
      		const precision cos_w0 = cos(omega);
      		const precision sin_w0 = sin(omega);

      		//Qualidade fixa (Butterworth 2ª ordem: ~0.707)
      		const precision quality = precision(1) / sqrt(precision(2));
      		const precision alpha   = sin_w0 / (precision(2) * quality);
      
      		//Coeficientes RBJ Low-pass
      		precision a0 = precision(1) + alpha;
  
      		precision b0 = (precision(1) - cos_w0) / (precision(2) * a0);
      		precision b1 = (precision(1) - cos_w0) / a0;
      		precision b2 = (precision(1) - cos_w0) / (precision(2) * a0);
      		precision a1 = (-precision(2) * cos_w0) / a0;
      		precision a2 = (precision(1) - alpha) / a0;
  
      		//Forma as seções SOS
      		for (size_t index = 0; index < this->m_num_sos; index++)
      		    this->m_sos_sections[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}
	};

	template<typename precision, size_t order>
	class RBJCookbookHighpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* Construtor do filtro RBJ Audio-EQ-Cookbook do tipo passa-alta
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte para passa-alta */
		RBJCookbookHighpass(const precision& sample_rate, const precision& first_cut) {
			//Atribui a frequência de amostragem e a frequência de corte
			this->m_sample_rate = sample_rate;
			this->m_first_cutoff = first_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//O filtro RBJ não utiliza polos ou zeros analógicos
			return true;
		}

		auto convertPrototype() -> bool override {
			//O filtro RBJ não utiliza conversão de protótipo
			return true;
		}

		auto createSOS() -> bool override {
      		//Frequência de corte digital pré-distorcida
      		const precision omega  = precision(2) * M_PI * this->m_first_cutoff / this->m_sample_rate;
      		const precision cos_w0 = cos(omega);
      		const precision sin_w0 = sin(omega);
  
      		//Qualidade fixa (Butterworth 2ª ordem: ~0.707)
      		const precision quality = precision(1) / sqrt(precision(2));
      		const precision alpha = sin_w0 / (precision(2) * quality);
  
      		//Coeficientes RBJ High-pass
      		precision a0 = precision(1) + alpha;
  
      		precision b0 = (precision(1) + cos_w0) / (precision(2) * a0);
      		precision b1 = -(precision(1) + cos_w0) / a0;
      		precision b2 = (precision(1) + cos_w0) / (precision(2) * a0);
      		precision a1 = (-precision(2) * cos_w0) / a0;
      		precision a2 = (precision(1) - alpha) / a0;
  
      		//Forma as seções SOS
      		for (size_t index = 0; index < this->m_num_sos; index++)
      		    this->m_sos_sections[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}
	};

	template<typename precision, size_t order>
	class RBJCookbookBandpass : public IIRGenericDesign<precision, 0, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* Construtor do filtro RBJ Audio-EQ-Cookbook do tipo passa-alta
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte inferior para passa-faixa
		 * @param precision& [Hz] Frequência de corte superior para passa-faixa */
		RBJCookbookBandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut) {
			//Atribui a frequência de amostragem e a frequência de corte
			this->m_sample_rate = sample_rate;
			this->m_first_cutoff  = first_cut;
			this->m_second_cutoff = second_cut;

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//O filtro RBJ não utiliza polos ou zeros analógicos
			return true;
		}

		auto convertPrototype() -> bool override {
			//O filtro RBJ não utiliza conversão de protótipo
			return true;
		}

		auto createSOS() -> bool override {
			//Calcula a largura de banda e a frequência de corte analógica
			const precision bandwidth = precision(2) * M_PI * (this->m_second_cutoff -  this->m_first_cutoff) / this->m_sample_rate;
			const precision omega 		= precision(2) * M_PI * sqrt(this->m_first_cutoff * this->m_second_cutoff) / this->m_sample_rate;

      		//ω para RBJ (digital) obtido a partir da pré-distorção analógica
      		const precision alpha  = sin(omega) * sinh(log(precision(2))/precision(2) * bandwidth / omega);
			const precision cos_w0 = cos(omega);

      		//coeficientes analíticos do RBJ para lowpass
			precision a0 = precision(1) + alpha;

      		precision b0 = alpha / a0;
      		precision b1 = 0;
      		precision b2 = -alpha / a0;
      		precision a1 = -precision(2) * cos_w0 / a0;
      		precision a2 = (precision(1) - alpha) / a0;

			//Forma os biquads a partir dos coeficientes de segunda ordem
			for (size_t index = 0; index < this->m_num_sos; index++)
				this->m_sos_sections[index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}
	};

}

#endif
