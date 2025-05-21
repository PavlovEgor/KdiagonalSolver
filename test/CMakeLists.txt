# cmake_minimum_required(VERSION 2.8)	 # Проверка версии CMake.
# 										# Если версия установленной программы
# 										# старее указаной, произойдёт аварийный выход.

# project(KdiagonalSolver)			# Название проекта

# set(SOURCE_EXE penta_test.cpp)		# Установка переменной со списком исходников для исполняемого файла

# set(SOURCE_LIB ../src/PentadiagonalSolver.cpp)			# Тоже самое, но для библиотеки

# add_library(libKdiagonalSOlver SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

# include_directories(../include)

# add_executable(penta_test ${SOURCE_EXE})	# Создает исполняемый файл с именем main

# set(MYCOMPILE_FLAGS "-lm")

# add_definitions(${MYCOMPILE_FLAGS})

# target_link_libraries(penta_test libKdiagonalSOlver)		# Линковка программы с библиотекой


cmake_minimum_required(VERSION 2.8)	 # Проверка версии CMake.
										# Если версия установленной программы
										# старее указаной, произойдёт аварийный выход.

project(KdiagonalSolver)			# Название проекта

set(SOURCE_EXE k_test.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB ../src/KdiagonalSolver.cpp)			# Тоже самое, но для библиотеки

add_library(libKdiagonalSOlver SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

include_directories(../include)

add_executable(k_test ${SOURCE_EXE})	# Создает исполняемый файл с именем main

set(MYCOMPILE_FLAGS "-lm")

add_definitions(${MYCOMPILE_FLAGS})

target_link_libraries(k_test libKdiagonalSOlver)		# Линковка программы с библиотекой