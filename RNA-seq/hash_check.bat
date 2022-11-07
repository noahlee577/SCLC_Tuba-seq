@echo off
for /R . %%f in (*.gz) do (
    echo | set/p="%%f - "
    certutil -hashfile "%%f" MD5 | findstr /V ":"
    cat "%%f.md5"
)