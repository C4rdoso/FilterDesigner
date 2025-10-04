#ifndef BUTTERWORTH_DESIGN_HPP
#define BUTTERWORTH_DESIGN_HPP

#include "filter_design.hpp"

namespace FilterDesigner {
	template<typename precision, size_t order>
	class ButterworthLowpass : public IIRGenericDesign<precision, order, 0, (order / 2) + (order % 2 == 1)> {
	public:
		/* Construtor do filtro Butterworth do tipo passa-baixa
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte para passa-baixa */
		ButterworthLowpass(const precision& sample_rate, const precision& first_cut) {
			//Atribui a frequência de amostragem
			this->m_sample_rate = sample_rate;

			//Pré-distorce as frequências para o domínio analógico
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);
			
			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calcula todos os polos e seus conjugados para o protótipo do filtro passa-baixa
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Armazena o polo calculado
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Armazena o conjugado logo na sequência
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//Se houver um polo real (em caso de ordem impar)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Determina a frequência de corte analógica
			const precision omega = this->m_first_cutoff;

			//O número de polos analógicos é igual à ordem do filtro
			const size_t num_analogue_poles = order;

			//Calcula os polos digitais, multiplicando-os pela frequência de corte
			for (size_t index = 0; index < num_analogue_poles; index++) {
				this->m_poles[index] *= omega;

				//Avalia se todos os polos estão no semiplano esquerdo
				if (this->m_poles[index].real() > 0) {
					//Caso contrário, interrompe o processo e gera uma excessão
					return false;
				}
			}

			//Calcula o ganho inicial do filtro
			precision initial_gain = pow(omega, num_analogue_poles);

			//Faz a conversão do plano, calcula os zeros e obtém o ganho após a conversão
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Faz uma correção do ganho
			this->m_overall_gain = initial_gain * (initial_gain / this->m_overall_gain);

			return true;
		}

		auto createSOS() -> bool override {
			//Gera as seções de segunda ordem a partir dos polos e zeros
			zpk2Biquads(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_sos_sections, this->m_num_sos
			);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}

	};

	template<typename precision, size_t order>
	class ButterworthHighpass : public IIRGenericDesign<precision, order, order, (order / 2) + (order % 2 == 1)> {
	public:
		/* Construtor do filtro Butterworth do tipo passa-alta
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte para passa-alta */
		ButterworthHighpass(const precision& sample_rate, const precision& first_cut) {
			//Atribui a frequência de amostragem
			this->m_sample_rate = sample_rate;

			//Pré-distorce as frequências para o domínio analógico
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);

			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calcula todos os polos e seus conjugados para o protótipo do filtro passa-baixa
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Armazena o polo calculado
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Armazena o conjugado logo na sequência
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//Se houver um polo real (em caso de ordem impar)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Determina a frequência de corte analógica
			const precision omega = this->m_first_cutoff;

			//O número de polos analógicos é igual à ordem do filtro
			const size_t num_analogue_poles = order;

			//Calcula os polos digitais, multiplicando-os pela frequência de corte
			for (size_t index = 0; index < num_analogue_poles; index++) {
				//Ignora polos nulos
				if (cabs(this->m_poles[index]) == 0) continue;

				//Converte o polo passa-baixa para passa-alta usando a frequência de corte
				this->m_poles[index] = Complex<precision>(omega) / this->m_poles[index];

				//Inicializa os zeros
				this->m_zeros[index] = Complex<precision>(0.0);

				//Avalia se todos os polos estão no semiplano esquerdo
				if (this->m_poles[index].real() > 0) {
					//Caso contrário, interrompe o processo e gera uma excessão
					return false;
				}
			}

			//Calcula o ganho inicial do filtro
			precision initial_gain = 1.0;

			//Faz a conversão do plano, calcula os zeros e obtém o ganho após a conversão
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Faz uma correção do ganho
			this->m_overall_gain = precision(1) / this->m_overall_gain;

			return true;
		}

		auto createSOS() -> bool override {
			//Gera as seções de segunda ordem a partir dos polos e zeros
			zpk2Biquads<precision>(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_sos_sections, this->m_num_sos
			);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}

	};

	template<typename precision, size_t order>
	class ButterworthBandpass : public IIRGenericDesign<precision, order * 2, order, order> {
	public:
		/* Construtor do filtro Butterworth do tipo passa-faixa
		 * @param precision& [Hz] Taxa de amostragem do sinal de entrada
		 * @param precision& [Hz] Frequência de corte inferior para passa-faixa
		 * @param precision& [Hz] Frequência de corte superior para passa-faixa */
		ButterworthBandpass(const precision& sample_rate, const precision& first_cut, const precision& second_cut) {
			//Atribui a frequência de amostragem
			this->m_sample_rate = sample_rate;
			
			//Pré-distorce as frequências para o domínio analógico
			this->m_first_cutoff = 2 * tan(M_PI * first_cut / sample_rate);
			this->m_second_cutoff = 2 * tan(M_PI * second_cut / sample_rate);
			
			this->design();
		}

	protected:
		auto createPrototype() -> bool override {
			//Calcula todos os polos e seus conjugados para o protótipo do filtro passa-baixa
			for (size_t k = 0; k < order / 2; k++) {
				double theta = (double)(2 * k + 1) * M_PI / (2 * order);
				double real = -sin(theta);
				double imag = cos(theta);

				//Armazena o polo calculado
				this->m_poles[k * 2] = Complex<precision>(real, imag);

				//Armazena o conjugado logo na sequência
				this->m_poles[(k * 2) + 1] = Complex<precision>(real, -imag);
			}

			//Se houver um polo real (em caso de ordem impar)
			if (order % 2 == 1) {
				this->m_poles[order - 1] = Complex<precision>(-1.0, 0.0);
			}

			return true;
		}

		auto convertPrototype() -> bool override {
			//Calcula a largura de banda e a frequência de corte analógica
			const precision bandwidth = this->m_second_cutoff - this->m_first_cutoff;
			const precision omega = sqrt(this->m_first_cutoff * this->m_second_cutoff);

			//O número de polos analógicos é igual à ordem do filtro
			const size_t num_analogue_poles = order;

			//Transforma o protótipo passa-baixa em passa-faixa
			//O número de polos no filtro passa-faixa é 2n
			for (size_t index = 0; index < num_analogue_poles; index++) {
				//Ignora polos nulos
				if (cabs(this->m_poles[index]) == 0) continue;

				//Caso seja um polo real (ultimo polo no filtro de ordem impar)
				if (this->m_poles[index].imag() == 0 && index == num_analogue_poles - 1 && (order % 2 == 1)) {
					Complex<precision> first_half  = precision(0.5) * this->m_poles[index] * bandwidth + Complex<precision>(0, omega);
					Complex<precision> second_half = precision(0.5) * this->m_poles[index] * bandwidth - Complex<precision>(0, omega);

					this->m_poles[index] = first_half;
					this->m_poles[index + num_analogue_poles] = second_half;

					//Avalia se todos os polos estão no semiplano esquerdo
					if (first_half.real() > 0 || second_half.real() > 0) {
						//Caso contrário, interrompe o processo e gera uma excessão
						return false;
					}
				}

				//Caso seja um polo complexo
				else {
					//Calcula o primeiro e o segundo termo
					Complex<precision> first  = precision(0.5) * this->m_poles[index] * bandwidth;
					Complex<precision> second = precision(0.5) * csqrt((bandwidth * bandwidth) * (this->m_poles[index] * this->m_poles[index]) - 4 * omega * omega);

					//A primeira metade deve ser o produto do primeiro com o segundo
					const auto& first_half_pole = this->m_poles[index] = first + second;

					//A segunda metade deve ser a diferença do primeiro com o segundo
					const auto& second_half_pole = this->m_poles[index + num_analogue_poles] = first - second;

					//Avalia se todos os polos estão no semiplano esquerdo
					if (first_half_pole.real() > 0 || second_half_pole.real() > 0) {
						//Caso contrário, interrompe o processo e gera uma excessão
						return false;
					}
				}
			}

			//Calcula o ganho inicial do filtro
			precision initial_gain = pow(bandwidth, num_analogue_poles);

			//Faz a conversão do plano, calcula os zeros e obtém o ganho após a conversão
			this->m_overall_gain = planeConversion(initial_gain, this->m_poles, this->m_num_poles, this->m_zeros, this->m_num_zeros);

			//Faz uma correção do ganho
			this->m_overall_gain = initial_gain * (initial_gain / this->m_overall_gain);

			return true;
		}

		auto createSOS() -> bool override {
			//Gera as seções de segunda ordem a partir dos polos e zeros
			zpk2Biquads<precision>(
				this->m_poles, this->m_num_poles,
				this->m_zeros, this->m_num_zeros,
				this->m_overall_gain,
				this->m_sos_sections, this->m_num_sos
			);

			//Caso tenha chegado até aqui, indica sucesso
			return true;
		}

	};
}
#endif
