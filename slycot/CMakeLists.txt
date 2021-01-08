# CMakeLists.txt file for Slycot
# use in conjunction with scikit-build
#
# RvP, 180710


set(FSOURCES

    src/AB01MD.f  src/AB01ND.f  src/AB01OD.f  src/AB04MD.f  src/AB05MD.f
    src/AB05ND.f  src/AB05OD.f  src/AB05PD.f  src/AB05QD.f  src/AB05RD.f
    src/AB05SD.f  src/AB07MD.f  src/AB07ND.f  src/AB08MD.f  src/AB08MZ.f
    src/AB08ND.f  src/AB08NX.f  src/AB08NZ.f  src/AB09AD.f  src/AB09AX.f
    src/AB09BD.f  src/AB09BX.f  src/AB09CD.f  src/AB09CX.f  src/AB09DD.f
    src/AB09ED.f  src/AB09FD.f  src/AB09GD.f  src/AB09HD.f  src/AB09HX.f
    src/AB09HY.f  src/AB09ID.f  src/AB09IX.f  src/AB09IY.f  src/AB09JD.f
    src/AB09JV.f  src/AB09JW.f  src/AB09JX.f  src/AB09KD.f  src/AB09KX.f
    src/AB09MD.f  src/AB09ND.f  src/AB13AD.f  src/AB13AX.f  src/AB13BD.f
    src/AB13CD.f  src/AB13DD.f  src/AB13DX.f  src/AB13ED.f  src/AB13FD.f
    src/AB13MD.f  src/AB8NXZ.f  src/AG07BD.f  src/AG08BD.f  src/AG08BY.f
    src/AG08BZ.f  src/AG8BYZ.f  src/BB01AD.f  src/BB02AD.f  src/BB03AD.f
    src/BB04AD.f  src/BD01AD.f  src/BD02AD.f  src/DE01OD.f  src/DE01PD.f
    src/DF01MD.f  src/DG01MD.f  src/DG01ND.f  src/DG01NY.f  src/DG01OD.f
    src/DK01MD.f  src/FB01QD.f  src/FB01RD.f  src/FB01SD.f  src/FB01TD.f
    src/FB01VD.f  src/FD01AD.f  src/IB01AD.f  src/IB01BD.f  src/IB01CD.f
    src/IB01MD.f  src/IB01MY.f  src/IB01ND.f  src/IB01OD.f  src/IB01OY.f
    src/IB01PD.f  src/IB01PX.f  src/IB01PY.f  src/IB01QD.f  src/IB01RD.f
    src/IB03AD.f  src/IB03BD.f  src/MA01AD.f  src/MA02AD.f  src/MA02BD.f
    src/MA02BZ.f  src/MA02CD.f  src/MA02CZ.f  src/MA02DD.f  src/MA02ED.f
    src/MA02FD.f  src/MA02GD.f  src/MA02HD.f  src/MA02ID.f  src/MA02JD.f
    src/MB01MD.f  src/MB01ND.f  src/MB01PD.f  src/MB01QD.f  src/MB01RD.f
    src/MB01RU.f  src/MB01RW.f  src/MB01RX.f  src/MB01RY.f  src/MB01SD.f
    src/MB01TD.f  src/MB01UD.f  src/MB01UW.f  src/MB01UX.f  src/MB01VD.f
    src/MB01WD.f  src/MB01XD.f  src/MB01XY.f  src/MB01YD.f  src/MB01ZD.f
    src/MB02CD.f  src/MB02CU.f  src/MB02CV.f  src/MB02CX.f  src/MB02CY.f
    src/MB02DD.f  src/MB02ED.f  src/MB02FD.f  src/MB02GD.f  src/MB02HD.f
    src/MB02ID.f  src/MB02JD.f  src/MB02JX.f  src/MB02KD.f  src/MB02MD.f
    src/MB02ND.f  src/MB02NY.f  src/MB02OD.f  src/MB02PD.f  src/MB02QD.f
    src/MB02QY.f  src/MB02RD.f  src/MB02RZ.f  src/MB02SD.f  src/MB02SZ.f
    src/MB02TD.f  src/MB02TZ.f  src/MB02UD.f  src/MB02UU.f  src/MB02UV.f
    src/MB02VD.f  src/MB02WD.f  src/MB02XD.f  src/MB02YD.f  src/MB03MD.f
    src/MB03MY.f  src/MB03ND.f  src/MB03NY.f  src/MB03OD.f  src/MB03OY.f
    src/MB03PD.f  src/MB03PY.f  src/MB03QD.f  src/MB03QX.f  src/MB03QY.f
    src/MB03RD.f  src/MB03RX.f  src/MB03RY.f  src/MB03SD.f  src/MB03TD.f
    src/MB03TS.f  src/MB03UD.f  src/MB03VD.f  src/MB03VY.f  src/MB03WA.f
    src/MB03WD.f  src/MB03WX.f  src/MB03XD.f  src/MB03XP.f  src/MB03XU.f
    src/MB03YA.f  src/MB03YD.f  src/MB03YT.f  src/MB03ZA.f  src/MB03ZD.f
    src/MB04DD.f  src/MB04DI.f  src/MB04DS.f  src/MB04DY.f  src/MB04GD.f
    src/MB04ID.f  src/MB04IY.f  src/MB04IZ.f  src/MB04JD.f  src/MB04KD.f
    src/MB04LD.f  src/MB04MD.f  src/MB04ND.f  src/MB04NY.f  src/MB04OD.f
    src/MB04OW.f  src/MB04OX.f  src/MB04OY.f  src/MB04PA.f  src/MB04PB.f
    src/MB04PU.f  src/MB04PY.f  src/MB04QB.f  src/MB04QC.f  src/MB04QF.f
    src/MB04QU.f  src/MB04TB.f  src/MB04TS.f  src/MB04TT.f  src/MB04TU.f
    src/MB04TV.f  src/MB04TW.f  src/MB04TX.f  src/MB04TY.f  src/MB04UD.f
    src/MB04VD.f  src/MB04VX.f  src/MB04WD.f  src/MB04WP.f  src/MB04WR.f
    src/MB04WU.f  src/MB04XD.f  src/MB04XY.f  src/MB04YD.f  src/MB04YW.f
    src/MB04ZD.f  src/MB05MD.f  src/MB05MY.f  src/MB05ND.f  src/MB05OD.f
    src/MB05OY.f  src/MB3OYZ.f  src/MB3PYZ.f  src/MC01MD.f  src/MC01ND.f
    src/MC01OD.f  src/MC01PD.f  src/MC01PY.f  src/MC01QD.f  src/MC01RD.f
    src/MC01SD.f  src/MC01SW.f  src/MC01SX.f  src/MC01SY.f  src/MC01TD.f
    src/MC01VD.f  src/MC01WD.f  src/MC03MD.f  src/MC03ND.f  src/MC03NX.f
    src/MC03NY.f  src/MD03AD.f  src/MD03BA.f  src/MD03BB.f  src/MD03BD.f
    src/MD03BF.f  src/MD03BX.f  src/MD03BY.f  src/NF01AD.f  src/NF01AY.f
    src/NF01BA.f  src/NF01BB.f  src/NF01BD.f  src/NF01BE.f  src/NF01BF.f
    src/NF01BP.f  src/NF01BQ.f  src/NF01BR.f  src/NF01BS.f  src/NF01BU.f
    src/NF01BV.f  src/NF01BW.f  src/NF01BX.f  src/NF01BY.f  src/SB01BD.f
    src/SB01BX.f  src/SB01BY.f  src/SB01DD.f  src/SB01FY.f  src/SB01MD.f
    src/SB02CX.f  src/SB02MD.f  src/SB02MR.f  src/SB02MS.f  src/SB02MT.f
    src/SB02MU.f  src/SB02MV.f  src/SB02MW.f  src/SB02ND.f  src/SB02OD.f
    src/SB02OU.f  src/SB02OV.f  src/SB02OW.f  src/SB02OX.f  src/SB02OY.f
    src/SB02PD.f  src/SB02QD.f  src/SB02RD.f  src/SB02RU.f  src/SB02SD.f
    src/SB03MD.f  src/SB03MU.f  src/SB03MV.f  src/SB03MW.f  src/SB03MX.f
    src/SB03MY.f  src/SB03OD.f  src/SB03OR.f  src/SB03OT.f  src/SB03OU.f
    src/SB03OV.f  src/SB03OY.f  src/SB03PD.f  src/SB03QD.f  src/SB03QX.f
    src/SB03QY.f  src/SB03RD.f  src/SB03SD.f  src/SB03SX.f  src/SB03SY.f
    src/SB03TD.f  src/SB03UD.f  src/SB04MD.f  src/SB04MR.f  src/SB04MU.f
    src/SB04MW.f  src/SB04MY.f  src/SB04ND.f  src/SB04NV.f  src/SB04NW.f
    src/SB04NX.f  src/SB04NY.f  src/SB04OD.f  src/SB04OW.f  src/SB04PD.f
    src/SB04PX.f  src/SB04PY.f  src/SB04QD.f  src/SB04QR.f  src/SB04QU.f
    src/SB04QY.f  src/SB04RD.f  src/SB04RV.f  src/SB04RW.f  src/SB04RX.f
    src/SB04RY.f  src/SB06ND.f  src/SB08CD.f  src/SB08DD.f  src/SB08ED.f
    src/SB08FD.f  src/SB08GD.f  src/SB08HD.f  src/SB08MD.f  src/SB08MY.f
    src/SB08ND.f  src/SB08NY.f  src/SB09MD.f  src/SB10AD.f  src/SB10DD.f
    src/SB10ED.f  src/SB10FD.f  src/SB10HD.f  src/SB10ID.f  src/SB10JD.f
    src/SB10KD.f  src/SB10LD.f  src/SB10MD.f  src/SB10PD.f  src/SB10QD.f
    src/SB10RD.f  src/SB10SD.f  src/SB10TD.f  src/SB10UD.f  src/SB10VD.f
    src/SB10WD.f  src/SB10YD.f  src/SB10ZD.f  src/SB10ZP.f  src/SB16AD.f
    src/SB16AY.f  src/SB16BD.f  src/SB16CD.f  src/SB16CY.f  src/SG02AD.f
    src/SG03AD.f  src/SG03AX.f  src/SG03AY.f  src/SG03BD.f  src/SG03BU.f
    src/SG03BV.f  src/SG03BW.f  src/SG03BX.f  src/SG03BY.f  src/TB01ID.f
    src/TB01IZ.f  src/TB01KD.f  src/TB01LD.f  src/TB01MD.f  src/TB01ND.f
    src/TB01PD.f  src/TB01TD.f  src/TB01TY.f  src/TB01UD.f  src/TB01VD.f
    src/TB01VY.f  src/TB01WD.f  src/TB01XD.f  src/TB01XZ.f  src/TB01YD.f
    src/TB01ZD.f  src/TB03AD.f  src/TB03AY.f  src/TB04AD.f  src/TB04AY.f
    src/TB04BD.f  src/TB04BV.f  src/TB04BW.f  src/TB04BX.f  src/TB04CD.f
    src/TB05AD.f  src/TC01OD.f  src/TC04AD.f  src/TC05AD.f  src/TD03AD.f
    src/TD03AY.f  src/TD04AD.f  src/TD05AD.f  src/TF01MD.f  src/TF01MX.f
    src/TF01MY.f  src/TF01ND.f  src/TF01OD.f  src/TF01PD.f  src/TF01QD.f
    src/TF01RD.f  src/TG01AD.f  src/TG01AZ.f  src/TG01BD.f  src/TG01CD.f
    src/TG01DD.f  src/TG01ED.f  src/TG01FD.f  src/TG01FZ.f  src/TG01HD.f
    src/TG01HX.f  src/TG01ID.f  src/TG01JD.f  src/TG01WD.f  src/UD01BD.f
    src/UD01CD.f  src/UD01DD.f  src/UD01MD.f  src/UD01MZ.f  src/UD01ND.f
    src/UE01MD.f

    src/delctg.f  src/select.f
    src/SLCT_DLATZM.f  src/SLCT_ZLATZM.f

    src/ftruefalse.f src/XERBLA.f
)

set(F2PYSOURCE src/_wrapper.pyf)
set(F2PYSOURCE_DEPS
  src/analysis.pyf src/math.pyf
  src/transform.pyf src/synthesis.pyf
  src/_helper.pyf)

configure_file(version.py.in version.py @ONLY)

set(PYSOURCE

  __init__.py  examples.py exceptions.py
  analysis.py math.py synthesis.py transform.py
  ${CMAKE_CURRENT_BINARY_DIR}/version.py)

set(SLYCOT_MODULE "_wrapper")

set(GENERATED_MODULE
  ${CMAKE_CURRENT_BINARY_DIR}/${SLYCOT_MODULE}${PYTHON_EXTENSION_MODULE_SUFFIX})


set(CMAKE_Fortran_FLAGS )

add_custom_target(wrapper ALL DEPENDS ${FSOURCES})
add_custom_command(
  OUTPUT _wrappermodule.c _wrapper-f2pywrappers.f
  COMMAND ${F2PY_EXECUTABLE} -m SLYCOT
  ${CMAKE_CURRENT_SOURCE_DIR}/${F2PYSOURCE}
  DEPENDS ${F2PYSOURCE_DEPS} ${F2PYSOURCE}
)

add_library(
  ${SLYCOT_MODULE} MODULE
  _wrappermodule.c
  ${F2PY_INCLUDE_DIR}/fortranobject.c
  _wrapper-f2pywrappers.f
  ${FSOURCES})

target_link_libraries(${SLYCOT_MODULE}
  ${LAPACK_LIBRARIES})

target_include_directories(
  ${SLYCOT_MODULE} PUBLIC
  ${F2PY_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  )

if (UNIX)
  if (APPLE)
    set_target_properties(${SLYCOT_MODULE} PROPERTIES
      LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${SLYCOT_MODULE} PROPERTIES
      LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

python_extension_module(${SLYCOT_MODULE})

install(TARGETS ${SLYCOT_MODULE} LIBRARY DESTINATION slycot)
install(FILES ${PYSOURCE} DESTINATION slycot)

add_subdirectory(tests)
