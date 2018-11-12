# Microsoft Developer Studio Project File - Name="Restud" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) QuickWin Application" 0x0107

CFG=Restud - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Restud.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Restud.mak" CFG="Restud - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Restud - Win32 Release" (based on "Win32 (x86) QuickWin Application")
!MESSAGE "Restud - Win32 Debug" (based on "Win32 (x86) QuickWin Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Restud - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /libs:qwin /nologo /warn:nofileopt
# ADD F90 /compile_only /libs:qwin /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"
# ADD LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"

!ELSEIF  "$(CFG)" == "Restud - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /libs:qwin /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /libs:qwin /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x407 /d "_DEBUG"
# ADD RSC /l 0x407 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /incremental:no /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Restud - Win32 Release"
# Name "Restud - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\Autval.f90
NODEP_F90_AUTVA=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\autvalsav.f90
NODEP_F90_AUTVAL=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\autvalsavtr.f90
NODEP_F90_AUTVALS=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Autvaltr.f90
NODEP_F90_AUTVALT=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BORROW.F90
NODEP_F90_BORRO=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BORROWTR.F90
NODEP_F90_BORROW=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\cm.f90
NODEP_F90_CM_F9=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Condist.f90
NODEP_F90_CONDI=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\DGRID.F90
NODEP_F90_DGRID=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\DGRID0.F90
NODEP_F90_DGRID0=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Dgrid2.f90
NODEP_F90_DGRID2=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\DISS.F90
NODEP_F90_DISS_=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Disstr.f90
NODEP_F90_DISST=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Enforce.f90
NODEP_F90_ENFOR=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\enforcetr.f90
NODEP_F90_ENFORC=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\exdem.f90
NODEP_F90_EXDEM=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FinSS.f90
NODEP_F90_FINSS=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Fputil.f
# End Source File
# Begin Source File

SOURCE=.\GINI2.F90
NODEP_F90_GINI2=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Guess.f90
NODEP_F90_GUESS=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Guess2.f90
NODEP_F90_GUESS2=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GUESSTR.F90
NODEP_F90_GUESST=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\IniSS.f90
NODEP_F90_INISS=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Main.f90
DEP_F90_MAIN_=\
	{$(INCLUDE)}"MSIMSL.mod"\
	
NODEP_F90_MAIN_=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Markov.f
NODEP_F90_MARKO=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Mctransnew.f90
DEP_F90_MCTRA=\
	{$(INCLUDE)}"MSIMSL.mod"\
	
NODEP_F90_MCTRA=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Newtonss.f90
NODEP_F90_NEWTO=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Params.f90
# End Source File
# Begin Source File

SOURCE=.\Residss.f90
NODEP_F90_RESID=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Resource.f90
NODEP_F90_RESOU=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Saving.f90
NODEP_F90_SAVIN=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Sysbin.f90
NODEP_F90_SYSBI=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Sysbintr.f90
NODEP_F90_SYSBIN=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Sysex.f90
NODEP_F90_SYSEX=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Sysextr.f90
NODEP_F90_SYSEXT=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\transition.f90
NODEP_F90_TRANS=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Valuecm.f90
NODEP_F90_VALUE=\
	".\Debug\params.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\WGINI.F90
NODEP_F90_WGINI=\
	".\Debug\params.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
