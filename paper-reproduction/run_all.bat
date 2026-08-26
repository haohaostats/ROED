@echo off
cd /d "%~dp0"
python run.py --stage all %*
if errorlevel 1 pause
