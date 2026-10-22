@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    py -3.12 -m venv .venv
    if errorlevel 1 (
        echo Python 3.12 is required. Install it, then run this launcher again.
        echo winget install --exact --id Python.Python.3.12 --source winget --scope user
        pause
        exit /b 1
    )
)
".venv\Scripts\python.exe" -m pip install --disable-pip-version-check -r requirements.lock
if errorlevel 1 (
    echo Dependency installation failed. Copy the error into our chat.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" -m streamlit run app.py --server.headless false --server.address 127.0.0.1
pause
