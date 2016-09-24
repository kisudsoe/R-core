@echo off
path C:\Program Files\R\R-3.3.1\bin; %path%
:path C:\Program Files\R\R-3.2.5\bin; %path%

: current path
pushd %~dp0

: set title
title correlation coefficient batch by Kim Seung-Soo

set csvcount=0
if exist *.csv (
  : make 'result' folder
  if not exist result md result

  for %%a in (*.csv) do (
    set /a csvcount=!scvcount!+1
    rscript corr_pval.cmd.r "%%a" "result\%%a"

    if errorlevel 1 (
      echo [LOG] Fiail to execution
      pause
    ) else (
      echo [LOG] Success to execution
    )
  )
)

if %csvcount%==0 (
  echo ==Error Message====================
  echo Not exist *.csv file in this folder
  echo Quit batch. Thank you.
  echo ===================================
  pause
)
endlocal
exit

:: 참고 ::
: 현재 실행중인 배치파일의 경로 얻기 https://sunhyeon.wordpress.com/2012/08/05/161/
: 배치파일 명령어 모음 http://snoopybox.co.kr/1404
: 배치파일 명령어 모음 http://mssp.tistory.com/entry/배치파일-명령어
: 배치파일 언어 규칙 http://jangpd007.tistory.com/163
: r 스크립트 cmd 실행 http://bhpark.tistory.com/61
