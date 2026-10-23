@echo off
if not defined DIA_CSC_DOTNET exit /b 2
if not defined DIA_CSC_DLL exit /b 2
"%DIA_CSC_DOTNET%" exec "%DIA_CSC_DLL%" %*
exit /b %errorlevel%
