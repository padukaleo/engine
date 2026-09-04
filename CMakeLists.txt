cmake_minimum_required(VERSION 3.22.1)
project("chem_engine")

# ── Flags optimasi untuk perangkat tua ───────────────
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CHEM_FLAGS
    -O2
    -fno-exceptions
    -fno-rtti
    -ffunction-sections
    -fdata-sections
)

# ── Source files ──────────────────────────────────────
add_library(
    chem_engine
    SHARED

    # Core
    Core/parser.cpp

    # Calc
    Core/calc/molar.cpp
    Core/calc/stoich.cpp
    Core/calc/solution.cpp
    Core/calc/gas.cpp
    Core/calc/thermo.cpp
    Core/calc/equilibrium.cpp
    Core/calc/acid_base.cpp
    Core/calc/electrochem.cpp
    Core/calc/colligative.cpp
    Core/calc/kinetics.cpp
    Core/calc/nuclear.cpp

    # Balancer
    Core/balancer/matrix.cpp
    Core/balancer/gaussian.cpp
    Core/balancer/balancer.cpp

    # JNI Bridge
    bridge/jni.cpp
)

# ── Include dirs ──────────────────────────────────────
target_include_directories(chem_engine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/Include/calc
    ${CMAKE_CURRENT_SOURCE_DIR}/Core
)

# ── Apply flags ───────────────────────────────────────
target_compile_options(chem_engine PRIVATE ${CHEM_FLAGS})

# ── Linker flags buang dead code ──────────────────────
target_link_options(chem_engine PRIVATE
    -Wl,--gc-sections
    -Wl,--strip-all
)

# ── Link library Android ──────────────────────────────
target_link_libraries(chem_engine
    android
    log
    m       # math.h (sqrt, log, pow, dll)
)