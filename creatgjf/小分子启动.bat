@echo off
chcp 65001 >nul
setlocal

REM 获取桌面路径
set "desktop=%USERPROFILE%\Desktop"

REM 主文件夹名称
set "main_folder_name=ros"

REM 子文件夹名称列表
set "subfolders=TS;INT;单点能量"

REM 创建主文件夹路径
set "main_folder_path=%desktop%\%main_folder_name%"
if not exist "%main_folder_path%" (
    mkdir "%main_folder_path%"
)

REM 创建子文件夹
setlocal enabledelayedexpansion
for %%f in (%subfolders%) do (
    set "subfolder=%%f"
    set "subfolder=!subfolder:;=!"
    set "subfolder_path=%main_folder_path%\!subfolder!"
    if not exist "!subfolder_path!" (
        mkdir "!subfolder_path!"
    )
)
endlocal

echo 主文件夹 '%main_folder_name%' 和子文件夹已成功创建在桌面上。

pause
