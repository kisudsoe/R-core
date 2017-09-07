# R commands

v1.0_140310
...
v1.2_170907

[TOC]

# Session 관리 함수

```r
ls() # 기억된 변수의 목록을 표시
rm(x) # 변수 x를 삭제함
rm(list=ls()) # 기억된 모든 변수를 삭제
search() # search path

save.image() # 현재 작업공간을 저장함(Default path: C:\Program Files\R\rw1071) "*.RData" 파일이 생성됨
save.image("folder path") # R은 디렉토리 구분을 "\\" 또는 "/"으로 함에 주의

list.files() # 작업공간내 파일 리스트 보여줌
```


# I. 데이터 입출력 함수

## I-1. 입력 관련

```r
scan() # 직접 입력을 받기 위한 함수. 엔터키 입력시 종료
scan(n=m) # m개의 값을 입력
scan(what="") # 입력값을 문자로 지정
scan(sep=",") # ,를 기준으로 데이터를 분리하여 입력 받음

getwd() # 현재 작업공간을 확인
setwd("C:\Users\PC\Documents\R") # 작업공간을 설정/변경

read.table("aa.txt") # 데이터 프레임으로 aa.txt를 읽어옴
read.table("aa.txt", header=T, row.ranmes=1) # 헤더와 행이름 지정
read.table("aa.txt", sep=",") # 데이터를 콤마(,)로 구분함
read.table("aa.txt", dec=",") # 숫자의 자리 구분을 콤마(,)로 함
read.table("aa.txt", na.strings=".") # 결측값을 마침표(.)로 구분

read.delim("aa.txt") # Tap deliminated 형식을 읽어옴
ream.delim2("aa.txt") # 탭(Tab)으로 구분. 콤마(,)는 숫자의 자리를 구분함
read.csv("aa.csv",header=T) # aa.csv 파일을 읽어옴
read.csv2("aa.csv") # 세미콜론(;)으로 구분함. 마침표(.)는 숫자의 자리를 구분

cat(x) # x변수의 내용을 화면에 출력
cat(x, file="aa.txt") # aa.txt에 벡터 x 출력
cat(mat, file="aa2.txt") # matrix 데이터도 벡터 형태로만 출력됨
```

## I-2. 결측값 처리

```r
is.na(x)# 결측값 여부에 대한 논리벡터
complete.cases(x1,x2,...) # 주어진 벡터들의 대응하는 원소들이 모두 결측값이 아닌 때 "TRUE"임.

na.rm=# 통계함수에서 사용. "TRUE"이면 결측값을 제외하고 "FALSE"이면 NA값을 보냄
na.last=# 함수 sort()에서 결측값이 있을 때 "TRUE"일 경우 가장 뒤에 두고, "FALSE"일 때 가장 앖에 두며, "NA"일 때 값을 제외함.
na.print=# 함수 summary()와 print.default()에서 결측값을 표시하는 문자 또는 문자열을 지정
```

## I-3. 외부 데이터 입력 관련

```r
library(foreign) # foreign 패키지를 이용

read.spss("file path") # SPSS 파일 읽기
read.ssd("file path") # SAS 파일 읽기
read.mtp("file path") # Minitab 파일 읽기
read.dta("file path") # Stata 파일 읽기
```

## I-3. 출력 관련

```r
write(mat, file="aa.txt") # 작업공간에 aa.txt 이름과 테이블형태로 데이터를 출력
write.table(mat, "aa.csv", sep=",", row.names=TRUE) # 작업공간에 aa.csv 파일을 생성
```


# II. 일반 자료 함수

## II-1. 자료형 확인 및 표 구성, 그룹화, 재구성

```r
mode(x) # 데이터 x의 저장타입(mode)을 확인함
str(x) # 데이터의 각 열을 구성하는 자료의 mode를 확인
class(x) # 데이터 x의 저장 형태를 확인함
summary(x) # x의 자료요약

length(x) # 벡터 x의 길이(length)를 나타내게 함
dim(x) # 데이터 프레임의 행렬 크기를 표시

head(x) # 데이터의 위쪽 6행을 표시
head(x,n=m) # 데이터 x의 위쪽 m개 행을 표시
tail(x,n=m) # 데이터 x의 아래쪽 m개 행을 표시

table(f1, ...) # 데이터를 표 형식으로 요약
tapply(x, f, mean) # 데이터의 평균에 대한 표

factor(x) # 벡터를 요인으로 변환
factor(x, levels)# 데이터에 대한 수준을 지정
factor(x, labels)# 계산된 각 level에 대한 이름을 지정
factor(x, exclude)# 제외할 값을 지정, 기본값은 "NA"임
levels(f) = names# 수준의 새로운 이름

cut(x, breaks) # 연속형 변수를 주어진 구간으로 그룹화
cut(x, labels)# 각 그룹에 대한 이름을 지정
cut(x, right)# 구간의 오른쪽 점을 포함할 지에 대한 인수. "FALSE"일 때 왼쪽 점 포함.
```

## II-2. 벡터 관련

### 1) 내장 벡터

```r
letters # 알파벳이 저장되어 있는 벡터
month.name # 12달의 영어이름이 저장되어 있는 벡터
```

### 2) 벡터 일반

```r
rnorm(n, mean=x, sd=y) # mean=x, sd=y를 따르는 n개의 난수 생성
rep(x,n) # x를 n번 반복하는 벡터를 생성
min(x) # 최소값
max(x) # 최대값
range(x) # 범위
sum(x) # 벡터 원소들의 합
prod(x) # 벡터 원소들의 곱
median(x) # 중간값
mean(x) # 편균값
var(x) # 분산
sd(x) # 표준편차
cor(x, y) # 상관계수
quantile(x)# 4분위수
```

### 3) 벡터 논리

```r
which(x == "A") # 벡터 x에서 "A"의 인덱스 위치를 반환
match("A", x) # 벡터 x에서 "A"의 인뎃스 숫자를 반환
intersect(x, y) # x와 y간의 교집합을 반환
setdiff(x, y) # x에 대해 y의 차집합을 반환
union(x, y) # x와 y의 합집합을 반환
sort(x) # 오름차순 정렬
sort(x, decreasing=T) # 내림차순 정렬
sort(x, na.last=T) # 결측치인 NA를 가장 마지막으로 두고 정렬
order(x)
rank(x)
```

### 4) 벡터 데이터 생성 및 변환

```r
numeric(25) # 25개의 0(zero) 생성
character(25)# 25개의 " "(문자열) 생성
logical(25)# 25개의 FLASE 생성
seq(-4, 4, 0.1)# -4에서 4까지 0.1간격의 숫자열 생성
1:10# same as seq(1,10,1)
c(5, 7, 9, 13, 1:5)# 연결: 5 7 9 13 1 2 3 4 5
rep(1, 10)# 1을 열 번 반복
gl(3, 2, 12)# 3 가지 수준에서 2개의 블록으로 반복하여 총 길이 12인 벡터 생성 (즉, 1 1 2 2 3 3 1 1 2 2 3 3)

as.numeric(x)# 숫자형으로 변환
as.character(x)# 문자열로 변환
as.logical(x)# 논리형으로 변환
factor(x)# 벡터 x의 요인을 생성 <- ??
```

## II-3. 행렬 관련

```r
matrix(x, nrow=n, ncol=m) # x값을 가지는 n행 m열 행렬을 생성
apply(x, 1, fn) # x행렬의 각 행에 fn 함수 적용
apply(x, 2, fn) # x행렬의 각 열에 fn 함수 적용
t(x) # x행렬의 행과 열을 바꿈
x * y # 대칭원소간 곱셈을 수행
x%*%y # 행렬의 곱셈을 수행
x + y # 대칭 원소간 덧셈을 수행
x - y # 대칭 원소간 뺄셈
x / y # 대칭 원소간 나눗셈

x[n, m] # n행 m열의 원소
x[, m] # m열의 모든 원소
```

## II-4. 배열 관련

```r
 array(x, dim=c(l,m,n)) # x값을 가지는 l행 m열의 n개 배열을 생성
```

## II-5. 데이터 프레임 관련

```r
data.frame(x, y, z) # x, y, z벡터를 가지는 데이터 프레임 생성
as.data.frame(x) # x를 데이터 프레임형으로 바꿈

edit(x) # 데이터 프레임 x 를 수정하는 gui 실행
rbind(y) # 행방향 합치기
cbind(y) # 열방향 합치기

names(x) # 변수명을 반환
colnames(x) # 열의 이름 반환/지정
rownames(x) # 행의 이름 반환/지정
rownames(data) = data[,1] # data의 first column 값을 rowname으로 지정
nrow(x) # 행의 개수
ncol(x) # 열의 개수
dim(x) # 행, 열의 차원(dimension)을 반환
mean(x) # 모든 열의 평균을 구함

x$a # 데이터 프레임의 a열(벡터)을 추출하여 반환

attach(x) # 데이터 프레임의 각 변수를 메모리에 로드
detach(x) # 로드된 변수를 해제

data[1] = NULL # data의 first column을 삭제
```


# III. 논리 및 연산자, 함수, 인덱싱

## III-1. 사칙연산

```r
+, -, *, / # 사칙연산
^ # 제곱
%/% # 정수 나눔
%% # 나눈 나머지
```

## III-2. 논리 및 관계 연산자

```r
 == # 같다
 != # 같지 않다
 <,> # 작다, 크다
 <=, >= # 작거나 같다, 크거나 같다
 & # 논리곱(logical AND)
 | # 논리합(logical OR)
 ! # 논리부정(logical NOT)

is.na(x) # 데이터 x에서 결측치 여부 확인
```

## III-3. 함수

```r
log(x)# 자연로그 함수
log10(x)# 상용로그 함수
exp(x)# 지수함수
sin(x)# sin 함수
cos(x)
tan(x)
asin(x)# arcsin 함수
acos(x)
atan(x)

min(x)# 벡터에서 최소값
min(x, y, ...) # 전체 벡터 중 최소값
max(x)# 벡터에서 최대값
range(x)# same as c(min(x), max(x))
pmin(x1, x2, ...) # 벡터들에서 대칭 원소 간의 최소값
pmax(x1, x2, ...) # 벡터들에서 대칭 원소 간의 최대값
length(x) # 벡터의 길이
sum(complete.cases(x)) # 벡터에서 결측값을 제외한 원소의 개수
```

## III-4. 인덱싱

```r
x[1] # 첫번째 원소
x[1:5] # 1-5번째 원소
x[y<=30] # 논리를 만족하는 원소
x[sex=="male"] # 요인 변수를 만족하는 원소
i=c(2,3,5,7,11); x[i] # 숫자형 변수를 만족하는 원소
l=(y<=30); x[l] # 논리 변수를 만족하는 원소

m[4,] # 4번째 행
m[,3] # 3번째 열
dfr[dfr$var<=30,] # 조건을 만족하는 부분 선택
subset(dfr$var<=30) # same as dfr[dfr$var<=30,]

grep(pattern, x) # pattern을 벡터 x에서 검색, 인덱스를 반환함
gsub("before", "after", x) # 벡터 x에서 "before"를 "after"로 치환함
```


# IV. 그래프 함수

## IV-1. Graphical parameters

```r
par(...) # 그래픽 출력 변수 설정. ()안에 아래의 명령어를 쓴다.
```
* adj # 글자 정렬
* bg # 배경색
* bty # 외곽선 형태 지정
* cex # 심볼, 글씨 크기 cex.axis, cex.lab, cex.main, cex.sub
* col # 색지정 col.axis, col.lab, col.main, col.sub
* font # 글씨체
* las # 글씨 방향
* lty # 선형태
* lwd # 선굵기
* mar # 그래프 여백 지정
* mfcol # 출력 배열 형식 c(nr, nc)
* mfrow # 출력 배열 형식 (행우선)
* pch # 심볼 형태 지정
* ps # 글씨 크기(point)
* pty # 출력 영역 "s": 사각, "m": 최대
* tck, tcl # 축 마크
* xaxs, yaxs # 축 형태
* xaxt, yaxt # 축 표시

## IV-2. Hihg-level 관련

```r
x11() # 비어있는 그림창을 새로 만듬
dev.set(dev.prev()) # 이전 그림창을 활성화
dev.set(dev.next()) # 다음 그림창을 활성화

plot(x) # matrix 또는 data.frame으로 plot을 그림
plot(x1, x2) # x1과 x2를 각 축으로 dot plot을 그림
plot(x, n, m) # n에서 m 사이 범위에서 x를 그림
plot(data1, data2, xlab="x축 label", ylab="y축 label", main="제목", col=column, las=0/1/2/3)

barplot(data)
boxplot(data)

hist(data)# 히스토그램
hist(data, main="제목", xlab="x축", las=0/1/2/3, col="color", prob=TRUE/FALSE)

plotDensity(data)
pie(data)
persp(data) # perspective plot
```

## IV-3. Low-level 관련

```r
points(data)

lines(data, col="color")
ablines(x=x1)

text(data)
axis(data)
```

## IV-4. Trellis (격자) 관련

```r
xyplot
bwplot(data1, data2, data3, layout=c(column, rows, pages), main="제목")
histogram(data1, data2, data3, layout=c(column, rows, pages), main="제목")

par(mfrow=c(1,2)) # 1행 2열로 plot 파티션 나누기
```


# V. 통계 함수

```r
mean(x)# 평균
sd(x)# 표준편차
var(x)# 분산
median(x)# 중앙값
quantile(x, p) # 벡터 x에서 (100 x p)%에 해당하는 값
cor(x, y)# 상관계수
```

## V-1. 확률분포

### 정규분포

```r
dnorm(x)# 밀도함수를 계산
pnorm(x)# 누적분포확률 P(X<=x)를 계산
qnorm(x)# P(X<=x)=p를 만족하는 x를 구함
rnorm(x)# 정규분포를 따르는 n개의 난수를 생성
```

### 여러가지 확률분포

```r
pnorm(x, mean, sd)# 정규분포
plnorm(x, mean, sd)# 로그-정규분포
pt(x, df)# 스튜던트 t-분포
pf(x, n1, n2)# F 분포
pchisq(x, df)# 카이제곱 분포
pbinorm(x, n, p)# 이항분포
ppois(x, lambda)# 포아송분포
punif(x, rate)# 균일분포
pexp(x, rate)# 지수분포
pgamma(x, shape, scale) # 감마 분포
pbeta(x, a, b)# 베타분포
```
