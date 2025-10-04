#Especifica a versão mínima do CMake requerida para a compilação
cmake_minimum_required(VERSION 3.12)

set(NOME_PROJETO "FilterDesigner")              #Defina o nome do projeto
set(VERSAO_PROJETO 0.0.1)					    #Defina a versão do projeto

#=================================================================================
#   PROPRIEDADES DE PROJETO
#=================================================================================
#Especifica o nome do projeto e a sua versão
project(${NOME_PROJETO} VERSION ${VERSAO_PROJETO})

#Versão da linguagem (padrão C++20)
set(CMAKE_CXX_STANDARD 20)

#Altera o diretório padrão de arquivos binários gerados pelo CMake
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

#Altera o diretório padrão de arquivos de biblioteca dinâmica gerados pelo CMake
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)

#Altera o diretório padrão de arquivos de biblioteca estática gerados pelo CMake
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)

#Diretiva especial para Microsoft Visual C++
if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

#=================================================================================
#   COMPILAÇÃO DO PROJETO
#=================================================================================
#Bibliotecas necessárias para a compilção 

#Procura e lista todos os arquivos de fonte .cpp dentro da pastar src
file(GLOB ARQUIVOS_FONTE CONFIGURE_DEPENDS "src/*.cpp")

#Cria a biblioteca estática passando os arquivos de fonte
add_executable(${NOME_PROJETO} ${ARQUIVOS_FONTE})

#Vincula na biblioteca o diretório de arquivos de inclusão
target_include_directories(${NOME_PROJETO} PUBLIC include)