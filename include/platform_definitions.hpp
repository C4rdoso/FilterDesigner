#ifndef PLATFORM_DEFINITIONS_HPP
#define PLATFORM_DEFINITIONS_HPP

#if defined(_WIN32) || defined(_WIN64)
#   define PLATFORM_WINDOWS     1
#elif defined(__linux__)
#   define PLATFORM_LINUX       1
#elif defined(__APPLE__)
#   define PLATFORM_DARWIN      1
#elif defined(ESP_PLATFORM)
#   define PLATFORM_ESPRESSIF   1
#elif defined(STM32F1xx) || defined(STM32F4xx) || defined(STM32F7xx) || defined(STM32H7xx) || defined(__STM32__)
#   define PLATFORM_STM32       1
#elif defined(__AVR__)
#   define PLATFORM_AVR         1
#else
#   error "Plataforma nao identificada. Adicione o suporte em platform_definitions.hpp"
#endif

//Identifica as plataformas que NÃO possuem suporte completo a biblioteca STL
#if defined(PLATFORM_AVR)
#   define HAS_STL_SUPPORT 0
#   warning "Bibliotecas STL nao suportadas nessa plataforma! Usando ferramentas substitutas."
#else
#   define HAS_STL_SUPPORT 1
#endif

#endif