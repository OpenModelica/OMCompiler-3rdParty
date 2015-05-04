:start

@echo # Configuration for Microsoft Build Engine > Makefile
@echo. >> Makefile
@echo. 
@set prefix=
@set omp=0
@set intel=0
@set fortran=0
@set mpi=0
@set saamg=0
@set longlong=0
@set longdouble=0
@set debug=0
@set usercflags=
@set userfflags=
@set userldflags=

:again

@if "%1" == "-h" goto usage
@if "%1" == "--help" goto usage
@if "%1" == "--prefix" goto setprefix
@if "%1" == "--enable-omp" goto setomp
@if "%1" == "--enable-mpi" goto setmpi
@if "%1" == "--enable-intel" goto setintel
@if "%1" == "--enable-fortran" goto setfortran
@if "%1" == "--enable-saamg" goto setsaamg
@if "%1" == "--enable-longlong" goto setlonglong
@if "%1" == "--enable-longdouble" goto setlongdouble
@if "%1" == "--enable-debug" goto setdebug
@if "%1" == "--cflags" goto setusercflags
@if "%1" == "--fflags" goto setuserfflags
@if "%1" == "--ldflags" goto setuserldflags
@if "%1" == "" goto genmakefile

:usage

@echo Usage: configure [options]
@echo.	
@echo Options:
@echo.	
@echo.	--prefix PREFIX		Install Lis in directory PREFIX
@echo.	--enable-omp		Build with OpenMP library
@echo.	--enable-mpi		Build with MPI library
@echo.	--enable-intel		Use Intel Compiler
@echo.	--enable-fortran	Build Fortran interface
@echo.	--enable-saamg		Build SA-AMG preconditioner
@echo.	--enable-longlong	Build with 64bit integer support
@echo.	--enable-longdouble	Build with long double support
@echo.	--enable-debug		Build with debug mode
@echo.	--cflags FLAG		Pass FLAG to C compiler
@echo.	--fflags FLAG		Pass FLAG to Fortran compiler
@echo.	--ldflags FLAG		Pass FLAG to linker
@echo.	
@goto end

:setprefix

@shift
@set prefix=%1
@shift
@goto again

:setomp

@set omp=1
@shift
@goto again

:setmpi

@set mpi=1
@shift
@goto again

:setintel

@set intel=1
@shift
@goto again

:setfortran

@set intel=1
@set fortran=1
@shift
@goto again

:setsaamg

@set intel=1
@set fortran=1
@set saamg=1
@shift
@goto again

:setlonglong

@set longlong=1
@shift
@goto again

:setlongdouble

@set longdouble=1
@shift
@goto again

:setdebug

@set debug=1
@shift
@goto again

:setusercflags

@shift
@set usercflags=%1
@shift
@goto again

:setuserfflags

@shift
@set userfflags=%1
@shift
@goto again

:setuserldflags

@shift
@set userldflags=%1
@shift
@goto again

:genmakefile

@echo # Installation directory >> Makefile
@if not "(%prefix%)" == "()" (
@echo PREFIX=%prefix% >> Makefile
) else (
@echo PREFIX=.. >> Makefile
)
@echo. >> Makefile

@if (%omp%) == (1) (
@echo # Build with OpenMP library >> Makefile
@echo omp=1 >> Makefile
@echo. >> Makefile
@echo.	Build with OpenMP library		= yes
) else (
@echo.	Build with OpenMP library		= no
)

@if (%mpi%) == (1) (
@echo # Build with MPI library >> Makefile
@echo mpi=1 >> Makefile
@echo. >> Makefile
@echo.	Build with MPI library 			= yes
) else (
@echo.	Build with MPI library 			= no
)

@if (%intel%) == (1) (
@echo # Use Intel Compiler >> Makefile
@echo intel=1 >> Makefile
@echo. >> Makefile
@echo.	Use Intel Compiler			= yes
) else (
@echo.	Use Intel Compiler			= no
)

@if (%fortran%) == (1) (
@echo # Build Fortran interface >> Makefile
@echo fortran=1 >> Makefile
@echo. >> Makefile
@echo.	Build Fortran interface			= yes
) else (
@echo.	Build Fortran interface			= no
)

@if (%saamg%) == (1) (
@echo # Build SA-AMG preconditioner >> Makefile
@echo saamg=1 >> Makefile
@echo. >> Makefile
@echo.	Build SA-AMG preconditioner		= yes
) else (
@echo.	Build SA-AMG preconditioner		= no
)

@if (%longlong%) == (1) (
@echo # Build with 64bit integer support >> Makefile
@echo longlong=1 >> Makefile
@echo. >> Makefile
@echo.	Build with 64bit integer support	= yes
) else (
@echo.	Build with 64bit integer support	= no
)

@if (%longdouble%) == (1) (
@echo # Build with long double support >> Makefile
@echo longdouble=1 >> Makefile
@echo. >> Makefile
@echo.	Build with long double support		= yes
) else (
@echo.	Build with long double support		= no
)

@if (%debug%) == (1) (
@echo # Enable Debugging >> Makefile
@echo debug=1 >> Makefile
@echo. >> Makefile
@echo.	Build with debug mode			= yes
) else (
@echo.	Build with debug mode			= no
)

@echo # User-defined C flags >> Makefile
@if not "(%usercflags%)" == "()" (
@echo USER_CFLAGS = %usercflags% >> Makefile
) else (
@echo USER_CFLAGS = >> Makefile
)
@echo. >> Makefile

@echo # User-defined Fortran flags >> Makefile
@if not "(%userfflags%)" == "()" (
@echo USER_FFLAGS = %userfflags% >> Makefile
) else (
@echo USER_FFLAGS = >> Makefile
)
@echo. >> Makefile

@echo # User-defined linker flags >> Makefile
@if not "(%userldflags%)" == "()" (
@echo USER_LDFLAGS = %userldflags% >> Makefile
) else (
@echo USER_LDFLAGS = >> Makefile
)
@echo. >> Makefile
) 

@echo. >> Makefile
@type Makefile.in >> Makefile
@echo.

:end
