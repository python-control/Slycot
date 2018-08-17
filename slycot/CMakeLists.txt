# CMakeLists.txt file for Slycot
# use in conjunction with scikit-build
#
# RvP, 180710


set(FSOURCES

  src/AB01MD.f src/MA02AD.f src/MB03YT.f src/NF01BW.f src/SB10KD.f
  src/AB01ND.f src/MA02BD.f src/MB03ZA.f src/NF01BX.f src/SB10LD.f
  src/AB01OD.f src/MA02BZ.f src/MB03ZD.f src/NF01BY.f src/SB10MD.f
  src/AB04MD.f src/MA02CD.f src/MB04DD.f src/SB01BD.f src/SB10PD.f
  src/AB05MD.f src/MA02CZ.f src/MB04DI.f src/SB01BX.f src/SB10QD.f
  src/AB05ND.f src/MA02DD.f src/MB04DS.f src/SB01BY.f src/SB10RD.f
  src/AB05OD.f src/MA02ED.f src/MB04DY.f src/SB01DD.f src/SB10SD.f
  src/AB05PD.f src/MA02FD.f src/MB04GD.f src/SB01FY.f src/SB10TD.f
  src/AB05QD.f src/MA02GD.f src/MB04ID.f src/SB01MD.f src/SB10UD.f
  src/AB05RD.f src/MA02HD.f src/MB04IY.f src/SB02CX.f src/SB10VD.f
  src/AB05SD.f src/MA02ID.f src/MB04IZ.f src/SB02MD.f src/SB10WD.f
  src/AB07MD.f src/MA02JD.f src/MB04JD.f src/SB02MR.f src/SB10YD.f
  src/AB07ND.f src/MB01MD.f src/MB04KD.f src/SB02MS.f src/SB10ZD.f
  src/AB08MD.f src/MB01ND.f src/MB04LD.f src/SB02MT.f src/SB10ZP.f
  src/AB08MZ.f src/MB01PD.f src/MB04MD.f src/SB02MU.f src/SB16AD.f
  src/AB08ND.f src/MB01QD.f src/MB04ND.f src/SB02MV.f src/SB16AY.f
  src/AB08NX.f src/MB01RD.f src/MB04NY.f src/SB02MW.f src/SB16BD.f
  src/AB08NZ.f src/MB01RU.f src/MB04OD.f src/SB02ND.f src/SB16CD.f
  src/AB09AD.f src/MB01RW.f src/MB04OW.f src/SB02OD.f src/SB16CY.f
  src/AB09AX.f src/MB01RX.f src/MB04OX.f src/SB02OU.f src/select.f
  src/AB09BD.f src/MB01RY.f src/MB04OY.f src/SB02OV.f src/SG02AD.f
  src/AB09BX.f src/MB01SD.f src/MB04PA.f src/SB02OW.f src/SG03AD.f
  src/AB09CD.f src/MB01TD.f src/MB04PB.f src/SB02OX.f src/SG03AX.f
  src/AB09CX.f src/MB01UD.f src/MB04PU.f src/SB02OY.f src/SG03AY.f
  src/AB09DD.f src/MB01UW.f src/MB04PY.f src/SB02PD.f src/SG03BD.f
  src/AB09ED.f src/MB01UX.f src/MB04QB.f src/SB02QD.f src/SG03BU.f
  src/AB09FD.f src/MB01VD.f src/MB04QC.f src/SB02RD.f src/SG03BV.f
  src/AB09GD.f src/MB01WD.f src/MB04QF.f src/SB02RU.f src/SG03BW.f
  src/AB09HD.f src/MB01XD.f src/MB04QU.f src/SB02SD.f src/SG03BX.f
  src/AB09HX.f src/MB01XY.f src/MB04TB.f src/SB03MD.f src/SG03BY.f
  src/AB09HY.f src/MB01YD.f src/MB04TS.f src/SB03MU.f
  src/SLCT_DLATZM.f src/AB09ID.f src/MB01ZD.f src/MB04TT.f
  src/SB03MV.f src/SLCT_ZLATZM.f src/AB09IX.f src/MB02CD.f
  src/MB04TU.f src/SB03MW.f src/TB01ID.f src/AB09IY.f src/MB02CU.f
  src/MB04TV.f src/SB03MX.f src/TB01IZ.f src/AB09JD.f src/MB02CV.f
  src/MB04TW.f src/SB03MY.f src/TB01KD.f src/AB09JV.f src/MB02CX.f
  src/MB04TX.f src/SB03OD.f src/TB01LD.f src/AB09JW.f src/MB02CY.f
  src/MB04TY.f src/SB03OR.f src/TB01MD.f src/AB09JX.f src/MB02DD.f
  src/MB04UD.f src/SB03OT.f src/TB01ND.f src/AB09KD.f src/MB02ED.f
  src/MB04VD.f src/SB03OU.f src/TB01PD.f src/AB09KX.f src/MB02FD.f
  src/MB04VX.f src/SB03OV.f src/TB01TD.f src/AB09MD.f src/MB02GD.f
  src/MB04WD.f src/SB03OY.f src/TB01TY.f src/AB09ND.f src/MB02HD.f
  src/MB04WP.f src/SB03PD.f src/TB01UD.f src/AB13AD.f src/MB02ID.f
  src/MB04WR.f src/SB03QD.f src/TB01VD.f src/AB13AX.f src/MB02JD.f
  src/MB04WU.f src/SB03QX.f src/TB01VY.f src/AB13BD.f src/MB02JX.f
  src/MB04XD.f src/SB03QY.f src/TB01WD.f src/AB13CD.f src/MB02KD.f
  src/MB04XY.f src/SB03RD.f src/TB01XD.f src/AB13DD.f src/MB02MD.f
  src/MB04YD.f src/SB03SD.f src/TB01XZ.f src/AB13DX.f src/MB02ND.f
  src/MB04YW.f src/SB03SX.f src/TB01YD.f src/AB13ED.f src/MB02NY.f
  src/MB04ZD.f src/SB03SY.f src/TB01ZD.f src/AB13FD.f src/MB02OD.f
  src/MB05MD.f src/SB03TD.f src/TB03AD.f src/AB13MD.f src/MB02PD.f
  src/MB05MY.f src/SB03UD.f src/TB03AY.f src/AB8NXZ.f src/MB02QD.f
  src/MB05ND.f src/SB04MD.f src/TB04AD.f src/AG07BD.f src/MB02QY.f
  src/MB05OD.f src/SB04MR.f src/TB04AY.f src/AG08BD.f src/MB02RD.f
  src/MB05OY.f src/SB04MU.f src/TB04BD.f src/AG08BY.f src/MB02RZ.f
  src/MB3OYZ.f src/SB04MW.f src/TB04BV.f src/AG08BZ.f src/MB02SD.f
  src/MB3PYZ.f src/SB04MY.f src/TB04BW.f src/AG8BYZ.f src/MB02SZ.f
  src/MC01MD.f src/SB04ND.f src/TB04BX.f src/BB01AD.f src/MB02TD.f
  src/MC01ND.f src/SB04NV.f src/TB04CD.f src/BB02AD.f src/MB02TZ.f
  src/MC01OD.f src/SB04NW.f src/TB05AD.f src/BB03AD.f src/MB02UD.f
  src/MC01PD.f src/SB04NX.f src/TC01OD.f src/BB04AD.f src/MB02UU.f
  src/MC01PY.f src/SB04NY.f src/TC04AD.f src/BD01AD.f src/MB02UV.f
  src/MC01QD.f src/SB04OD.f src/TC05AD.f src/BD02AD.f src/MB02VD.f
  src/MC01RD.f src/SB04OW.f src/TD03AD.f src/DE01OD.f src/MB02WD.f
  src/MC01SD.f src/SB04PD.f src/TD03AY.f src/DE01PD.f src/MB02XD.f
  src/MC01SW.f src/SB04PX.f src/TD04AD.f src/delctg.f src/MB02YD.f
  src/MC01SX.f src/SB04PY.f src/TD05AD.f src/DF01MD.f src/MB03MD.f
  src/MC01SY.f src/SB04QD.f src/TF01MD.f src/DG01MD.f src/MB03MY.f
  src/MC01TD.f src/SB04QR.f src/TF01MX.f src/DG01ND.f src/MB03ND.f
  src/MC01VD.f src/SB04QU.f src/TF01MY.f src/DG01NY.f src/MB03NY.f
  src/MC01WD.f src/SB04QY.f src/TF01ND.f src/DG01OD.f src/MB03OD.f
  src/MC03MD.f src/SB04RD.f src/TF01OD.f src/DK01MD.f src/MB03OY.f
  src/MC03ND.f src/SB04RV.f src/TF01PD.f src/FB01QD.f src/MB03PD.f
  src/MC03NX.f src/SB04RW.f src/TF01QD.f src/FB01RD.f src/MB03PY.f
  src/MC03NY.f src/SB04RX.f src/TF01RD.f src/FB01SD.f src/MB03QD.f
  src/MD03AD.f src/SB04RY.f src/TG01AD.f src/FB01TD.f src/MB03QX.f
  src/MD03BA.f src/SB06ND.f src/TG01AZ.f src/FB01VD.f src/MB03QY.f
  src/MD03BB.f src/SB08CD.f src/TG01BD.f src/FD01AD.f src/MB03RD.f
  src/MD03BD.f src/SB08DD.f src/TG01CD.f src/IB01AD.f src/MB03RX.f
  src/MD03BF.f src/SB08ED.f src/TG01DD.f src/IB01BD.f src/MB03RY.f
  src/MD03BX.f src/SB08FD.f src/TG01ED.f src/IB01CD.f src/MB03SD.f
  src/MD03BY.f src/SB08GD.f src/TG01FD.f src/IB01MD.f src/MB03TD.f
  src/NF01AD.f src/SB08HD.f src/TG01FZ.f src/IB01MY.f src/MB03TS.f
  src/NF01AY.f src/SB08MD.f src/TG01HD.f src/IB01ND.f src/MB03UD.f
  src/NF01BA.f src/SB08MY.f src/TG01HX.f src/IB01OD.f src/MB03VD.f
  src/NF01BB.f src/SB08ND.f src/TG01ID.f src/IB01OY.f src/MB03VY.f
  src/NF01BD.f src/SB08NY.f src/TG01JD.f src/IB01PD.f src/MB03WA.f
  src/NF01BE.f src/SB09MD.f src/TG01WD.f src/IB01PX.f src/MB03WD.f
  src/NF01BF.f src/SB10AD.f src/UD01BD.f src/IB01PY.f src/MB03WX.f
  src/NF01BP.f src/SB10DD.f src/UD01CD.f src/IB01QD.f src/MB03XD.f
  src/NF01BQ.f src/SB10ED.f src/UD01DD.f src/IB01RD.f src/MB03XP.f
  src/NF01BR.f src/SB10FD.f src/UD01MD.f src/IB03AD.f src/MB03XU.f
  src/NF01BS.f src/SB10HD.f src/UD01MZ.f src/IB03BD.f src/MB03YA.f
  src/NF01BU.f src/SB10ID.f src/UD01ND.f src/MA01AD.f src/MB03YD.f
  src/NF01BV.f src/SB10JD.f src/UE01MD.f)

set(F2PYSOURCE src/_wrapper.pyf)

configure_file(version.py.in version.py @ONLY)

set(PYSOURCE

  __init__.py analysis.py examples.py math.py synthesis.py
  transform.py ${CMAKE_CURRENT_BINARY_DIR}/version.py)

set(SLYCOT_MODULE "_wrapper")
set(GENERATED_MODULE
  ${CMAKE_CURRENT_BINARY_DIR}/${SLYCOT_MODULE}${PYTHON_EXTENSION_MODULE_SUFFIX})


set(CMAKE_Fortran_FLAGS )

add_custom_target(wrapper ALL DEPENDS ${FSOURCES})
add_custom_command(
  OUTPUT SLYCOTmodule.c _wrappermodule.c _wrapper-f2pywrappers.f
  COMMAND ${F2PY_EXECUTABLE} -m SLYCOT
  ${CMAKE_CURRENT_SOURCE_DIR}/${F2PYSOURCE}
  )

add_library(
  ${SLYCOT_MODULE} SHARED
  SLYCOTmodule.c _wrappermodule.c _wrapper-f2pywrappers.f
  "${PYTHON_SITE}/numpy/f2py/src/fortranobject.c"
  ${FSOURCES})

set(CMAKE_SHARED_LIBRARY_PREFIX "")
if (WIN32)
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".pyd")
endif()
set_target_properties(${SLYCOT_MODULE} PROPERTIES
  OUTPUT_NAME "_wrapper")
if (WIN32)
  target_link_libraries(${SLYCOT_MODULE} PUBLIC
    ${PYTHON_LIBRARIES} ${LAPACK_LIBRARIES}) 
endif()

if (UNIX)
  target_link_libraries(${SLYCOT_MODULE} PUBLIC
    ${LAPACK_LIBRARIES})

  if (APPLE)
    set_target_properties(${SLYCOT_MODULE} PROPERTIES
      LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
    string(REGEX REPLACE "^([0-9]+)\.([0-9]+)\.[0-9]+$" "\\1\\2"
      PYMAJORMINOR ${PYTHON_VERSION_STRING})
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".cpython-${PYMAJORMINOR}m-darwin.so")
    message(STATUS "binary module suffix ${CMAKE_SHARED_LIBRARY_SUFFIX}")
  else()
    set_target_properties(${SLYCOT_MODULE} PROPERTIES
      LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

target_include_directories(
  ${SLYCOT_MODULE} PUBLIC
  ${PYTHON_SITE}/numpy/core/include
  ${PYTHON_SITE}/numpy/f2py/src
  ${PYTHON_INCLUDE_DIRS}
  )

install(TARGETS ${SLYCOT_MODULE} DESTINATION slycot)
install(FILES ${PYSOURCE} DESTINATION slycot)

add_subdirectory(tests)
