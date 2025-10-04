#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

//Inclui as definições padrões por plataforma
#include "platform_definitions.hpp"

#if HAS_STL_SUPPORT
#	include <cmath>		//Ferramentas matemáticas básicas da STL
#	include <cstdint>  	//Tipos inteiros com tamanho fixo
#	include <complex>	//Suporte a números complexos da STL
#endif

namespace FilterDesigner {

#if !HAS_STL_SUPPORT
	template <typename precision>
	struct Complex {
		//Parte real e imaginária do número complexo
		precision m_real;
		precision m_imag;

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

		//Operações Complexo <-> Real
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

		//Atribuições compostas (+=, -=, *=, /=) com Complex
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

		//Atribuições compostas com real
		auto operator+=(precision value) -> Complex& { m_real += value; return *this; }
		auto operator-=(precision value) -> Complex& { m_real -= value; return *this; }
		auto operator*=(precision value) -> Complex& { m_real *= value; m_imag *= value; return *this; }
		auto operator/=(precision value) -> Complex& { m_real /= value; m_imag /= value; return *this; }
		
		/* @return O valor real do número complexo */
		auto real() const -> precision { return m_real; }
		/* @return O valor imaginário do número complexo */
		auto imag() const -> precision { return m_imag; }
	};


	// --- Operadores globais Real <-> Complex (quando real está à esquerda) ---
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

	/* Calcula a magnetude ou norma do número complexo 
	 * @param Complex<precision>& Número complexo de entrada 
	 * @return A magnetude ou norma do número complexo */
	template<typename precision>
	auto cabs(const Complex<precision>& target) -> precision { 
		return sqrt(target.m_real * target.m_real + target.m_imag * target.m_imag); 
	}

	/* Calcula a raiz quadrada de um número complexo
	 * @param Complex<precision>& Número complexo de entrada 
	 * @return O número complexo resultante */
	template<typename precision>
	auto csqrt(const Complex<precision>& value) -> Complex<precision> {
		precision magnitude = sqrt(sqrt(value.m_real * value.m_real + value.m_imag * value.m_imag));
		precision angle = atan2(value.m_imag, value.m_real) / precision(2);
		return Complex<precision>(magnitude * cos(angle), magnitude * sin(angle));
	}

#else
	//Usa a implementação de números complexos da STL
	template<typename precision>
	using Complex = std::complex<precision>;

	//Apenas usa a função max implementada pela stl
	using std::max;

	/* Calcula a magnetude ou norma do número complexo 
	 * @param Complex<precision>& Número complexo de entrada 
	 * @return A magnetude ou norma do número complexo */
	template<typename precision>
	auto cabs(const Complex<precision>& target) -> precision { 
		return std::abs(target);
	}

	/* Calcula a raiz quadrada de um número complexo
	 * @param Complex<precision>& Número complexo de entrada 
	 * @return O número complexo resultante */
	template<typename precision>
	auto csqrt(const Complex<precision>& value) -> Complex<precision> {
		precision magnitude = std::sqrt(std::sqrt(value.real() * value.real() + value.imag() * value.imag()));
		precision angle = atan2(value.imag(), value.real()) / precision(2);
		return Complex<precision>(magnitude * cos(angle), magnitude * sin(angle));
	}

#endif

	template<typename precision>
	struct BiquadSection {
		//Coeficientes do filtro biquad de segunda ordem
		precision m_b0, m_b1, m_b2, m_a1, m_a2;

		//Controle de processo
		precision m_z1 = 0.0, m_z2 = 0.0;

		/* Construtor padrão da sessão de segunda ordem */
		BiquadSection()
		: m_b0(0.0), m_b1(0.0), m_b2(0.0), m_a1(0.0), m_a2(0.0) { }

		/* Cria uma sessão de segunda ordem na forma direta 2
		 * @param double& b Numeradores da sessão
		 * @param double& a Denominadores da sessão */
		BiquadSection(const precision& b0, const precision& b1, const precision& b2, const precision& a0, const precision& a1, const precision& a2)
		: m_b0(b0 / a0), m_b1(b1 / a0), m_b2(b2 / a0), m_a1((-a1) / a0), m_a2((-a2) / a0) { }

		/* Realiza o processamento de um sinal de entrada
		 * @param double& x Sinal de entrada
		 * @return O Sinal de entrada processado pela sessão */
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
		//O ganho após a conversão deve ser um ajuste do ganho inicial
		double blt_gain = initial_gain;

		//Ajusta o ganho e aplica a transformação bilinear dos zeros
		for (uint32_t i = 0; i < num_zeros; i++) {
			blt_gain /= bilinearTransform<precision>(zeros[i]);
		}

		//Ajusta o ganho e aplica a transformação bilinear dos polos
		for (uint32_t i = 0; i < num_poles; i++) {
			blt_gain *= bilinearTransform<precision>(poles[i]);
		}

		//Retorna o ganho após a transformação
		return blt_gain;
	}

	template<typename precision>
	inline auto zpk2Biquads(Complex<precision>* poles, const size_t num_poles, Complex<precision>* zeros, const size_t num_zeros, precision overall_gain, BiquadSection<precision>* biquads, const size_t num_biquads) -> size_t {
		//Ordem do filtro
		const size_t filter_order = max(num_zeros, num_poles);

		//Número de seções necessárias
		const size_t num_sections = (filter_order + 1) / 2;

		//Verifica se temos espaço suficiente no array
		if (num_sections > num_biquads) return false;

		size_t last_section_index = 0;

		//Processa pares
		for (size_t i = 0; i + 1 < filter_order; i += 2, last_section_index++) {
			const auto first_zero = i < num_zeros ? zeros[i] :			Complex<precision>(-1, 0);
			const auto second_zero = i + 1 < num_zeros ? zeros[i + 1] : Complex<precision>(-1, 0);
			const auto first_pole = i < num_poles ? poles[i] :			Complex<precision>(0, 0);
			const auto second_pole = i + 1 < num_poles ? poles[i + 1] : Complex<precision>(0, 0);

			//Numerador
			precision b0 = 1.0;
			precision b1 = -(first_zero + second_zero).real();
			precision b2 = (first_zero * second_zero).real();

			//Denominador
			precision a0 = 1.0;
			precision a1 = -(first_pole + second_pole).real();
			precision a2 = (first_pole * second_pole).real();

			//Aplica o ganho global apenas na primeira seção
			if (last_section_index == 0) {
				b0 *= overall_gain;
				b1 *= overall_gain;
				b2 *= overall_gain;
			}

			biquads[last_section_index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);
		}

		//Caso a ordem do filtro seja ímpar, o último polo e zero são tratados separadamente
		if (filter_order % 2 == 1) {
			const auto last_zero = (filter_order - 1 < num_zeros) ? zeros[filter_order - 1] : Complex<precision>(-1, 0);
			const auto last_pole = (filter_order - 1 < num_poles) ? poles[filter_order - 1] : Complex<precision>(0, 0);

			precision b0 = 1.0;
			precision b1 = -last_zero.real();
			precision b2 = 0.0;

			precision a0 = 1.0;
			precision a1 = -last_pole.real();
			precision a2 = 0.0;

			//Aplica ganho na primeira seção (se ela for justamente essa)
			if (last_section_index == 0) {
				b0 *= overall_gain;
				b1 *= overall_gain;
				b2 *= overall_gain;
			}

			biquads[last_section_index] = BiquadSection<precision>(b0, b1, b2, a0, -a1, -a2);
			last_section_index++;
		}

		//Indica para o usuário quantas sessões foram realmente geradas
		return last_section_index; 
	}

}

#endif
