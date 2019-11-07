# Install Jupyter Lab for R

`180921 FRI` 블로그에서 옮김

- R is already installed.

## 1. Upgrade R version

- In CMD (authorized), execute R.

```R
install.packages('stringr')
if(!require(installr)) { install.packages('installr'); require(installr) }
check.for.updates.R() # tells you if there is a new version of R or not.
```

- If there is no error message, execute following code line.

```R
updateR(F,T,T,F,T,F,T) # install,move,update.package,quit R
```

- 시스템 환경 변수 설정
- 윈도우 찾기: '시스템 환경 변수 편집' 
- Path에 등록된 옛날 R 주소를 새 버전의 주소로 바꾼다(주소 중간 버전 정보 수정).
- Or in CMD (authorized), execute following code

```cmd
set path=%path%;C:\Progarm Files\R\R-3.4.1\bin
```



## 2. Install Jupyter

- In CMD (authorized), execute following code for python version check

```cmd
python-V
```



- Execute codes for Jupyter and JupyterLab install
- JupyterLab requires Jupyter notebook version 4.3 or later.

```cmd
pip3 install --upgrade pip
pip3 install jupyter
pip3 install jupyterlab
```

```cmd
jupyter notebook --version
jupyter lab --version
```

- jupyter notebook: `5.6.0`, jupyter lab: `0.33.11`

- Project Jupyter homepage: <http://jupyter.org/install.html>



## 3. Install R kernel for Jupyter

- In R, install IRkernel. Be sure all JupyterLab was closed.

```R
install.packages(c('repr', 'IRdisplay', 'crayon', 'pbdZMQ', 'devtools','stringi','Rcpp'))
devtools::install_github('IRkernel/IRkernel')

IRkernel::installspec(name='ir',displayname='R') # 기존 R 커널을 대체함
IRkernel::installspec(user = F) # system-wide user available
```

- Per default IRkernel::installspec() will install a kernel with the name 'ir' and a display name of "R". (180921 current version is 0.7. Check the latest version by bellow link.)
  <https://irkernel.github.io/docs/IRkernel/>
- In Jupyter R session, check the session information by bellowing code.

```R
sessionInfo()
```

> ```
> R version 3.5.1 (2018-07-02)
> Platform: x86_64-w64-mingw32/x64 (64-bit)
> Running under: Windows 10 x64 (build 17134)
> 
> Matrix products: default
> 
> locale:
> [1] LC_COLLATE=Korean_Korea.949  LC_CTYPE=Korean_Korea.949   
> [3] LC_MONETARY=Korean_Korea.949 LC_NUMERIC=C                
> [5] LC_TIME=Korean_Korea.949    
> 
> attached base packages:
> [1] stats     graphics  grDevices utils     datasets  methods   base     
> 
> loaded via a namespace (and not attached):
>  [1] compiler_3.5.1       magrittr_1.5         IRdisplay_0.5.0     
>  [4] pbdZMQ_0.3-3         tools_3.5.1          htmltools_0.3.6     
>  [7] base64enc_0.1-3      crayon_1.3.4         Rcpp_0.12.18        
> [10] uuid_0.1-2           stringi_1.1.7        IRkernel_0.8.12.9000
> [13] jsonlite_1.5         stringr_1.3.1        digest_0.6.17       
> [16] repr_0.15.0          evaluate_0.11 
> ```