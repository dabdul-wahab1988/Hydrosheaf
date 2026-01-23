@echo off
setlocal

echo ==========================================
echo Hydrosheaf PEST++ Custom Builder
echo ==========================================

rem 1. Setup Environment
echo [1/4] Setting up Visual Studio 2022 environment...
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
if %errorlevel% neq 0 (
    echo Error: Failed to setup VS environment.
    exit /b 1
)

set CMAKE_EXE="C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
set NINJA_EXE="C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe"

rem 2. Prepare Build Directory
echo [2/4] Resuming build...
rem if exist pestpp\build rmdir /s /q pestpp\build
rem mkdir pestpp\build
cd pestpp\build

rem 3. Configure and Build
echo [3/4] Building with Ninja...
rem %CMAKE_EXE% -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=%NINJA_EXE% -DCMAKE_CXX_COMPILER=cl -DCMAKE_C_COMPILER=cl ..

%NINJA_EXE%
if %errorlevel% neq 0 (
    echo Error: CMake configuration failed.
    exit /b 1
)

%NINJA_EXE%
if %errorlevel% neq 0 (
    echo Error: Build failed.
    exit /b 1
)

rem 4. Install
echo [4/4] Installing binaries to bin/...
cd ..\..
if not exist bin mkdir bin

echo Copying executables...
for /r pestpp\build %%f in (*.exe) do (
    echo Found %%f
    copy "%%f" bin\ >nul
)

echo.
echo Build Complete!
echo Binaries are located in %CD%\bin
dir bin\*.exe
