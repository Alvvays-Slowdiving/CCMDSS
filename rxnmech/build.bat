@echo off
REM --- 이 스크립트는 rxnmech 프로젝트를 빌드하고 실행합니다. ---
REM 사용법:
REM   build.bat                    (기본값: 릴리즈(-O3) 모드로 빌드 및 실행)
REM   build.bat release              (릴리즈(-O3) 모드로 빌드 및 실행)
REM   build.bat debug                (디버그 모드로 빌드 및 실행)
REM   build.bat debug 50 40 30       (디버그 모드로 빌드 후, 인자를 전달하여 실행)

REM --- 변수 설정 ---
SET BUILD_DIR=build
SET CMAKE_GENERATOR="MinGW Makefiles"

REM 기본 빌드 타입을 Release로 설정
SET BUILD_TYPE=Release

REM --- 첫 번째 인자로 빌드 타입 결정 ---
IF /I "%1"=="debug" (
    SET BUILD_TYPE=Debug
    SHIFT /1
) ELSE IF /I "%1"=="release" (
    REM 'release'는 기본값이므로, 이 인자는 그냥 무시하고 넘어감
    SHIFT /1
)

echo =================================================
echo.            BUILD MODE: %BUILD_TYPE%
echo =================================================

REM --- 빌드 디렉토리 확인 및 생성 ---
IF NOT EXIST %BUILD_DIR% (
    echo Creating build directory: %BUILD_DIR%
    mkdir %BUILD_DIR%
)

REM --- CMake 실행 ---
echo.
echo Running CMake to configure the project...
cd %BUILD_DIR%
cmake .. -G %CMAKE_GENERATOR% -D CMAKE_BUILD_TYPE=%BUILD_TYPE%

REM CMake 실행 결과 확인
IF %ERRORLEVEL% NEQ 0 (
    echo.
    echo CMake configuration failed. Aborting.
    cd ..
    goto:eof
)

REM --- 컴파일 실행 (mingw32-make) ---
echo.
echo Compiling the project with mingw32-make...
mingw32-make -j%NUMBER_OF_PROCESSORS% --no-print-directory

REM 컴파일 결과 확인
IF %ERRORLEVEL% NEQ 0 (
    echo.
    echo Compilation failed. Aborting.
    cd ..
    goto:eof
)

REM --- 실행 ---
echo.
echo Build successful.
cd ..

REM build.bat에 전달된 나머지 모든 인자를 %* 로 받아서 실행 파일에 그대로 전달
echo Running the application: build\rxnmech_app.exe %*
echo --------------------------------------------------

REM --- 여기가 수정된 부분 ---
REM %BUILD_DIR% 변수 대신 직접 경로를 지정합니다.
build\rxnmech_app.exe %*

echo --------------------------------------------------
echo.
echo Application finished.
pause