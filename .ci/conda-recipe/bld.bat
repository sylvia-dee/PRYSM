xcopy /e "%RECIPE_DIR%\.." "%SRC_DIR%"
"%PYTHON%" setup.py install --single-version-externally-managed --record record.txt
if errorlevel 1 exit 1