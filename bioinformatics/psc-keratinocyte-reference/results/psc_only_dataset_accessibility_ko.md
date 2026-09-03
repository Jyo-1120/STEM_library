# PSC→keratinocyte 공개 데이터 접근성 재점검

점검일: 2026-09-03

## 결론

PSC 유래 데이터가 아예 없는 것은 아니다. processed/count matrix로 바로 분석 가능한 핵심 longitudinal 자료가 두 세트(GSE120107, GSE122383), 반복은 없지만 초기 시점이 촘촘한 보조 자료가 한 세트(GSE144241) 있다. scRNA-seq와 후기 endpoint 자료도 존재한다. 부족한 것은 **우리 2D RA/BMP4-EGF-CaCl2 프로토콜과 같은 조건에서 D0-D35를 촘촘하게 반복 측정한 단일 데이터셋**이다.

따라서 날짜 모델과 생물학적 marker 참고자료를 분리한다.

- 날짜 모델/검증: PSC 유래 자료만 사용
- primary keratinocyte 자료: terminal maturation marker의 생물학적 타당성 확인에만 사용
- 서로 다른 PSC 프로토콜: 원시 발현값을 직접 합치지 않고, study 내부 표준화 후 단계·방향·순서를 비교

## 추천 순위

| 우선순위 | 데이터셋 | PSC 유래 | 시점/반복 | 바로 쓸 수 있는 자료 | 원자료 접근성과 용량 | 권장 역할 |
|---|---|---:|---|---|---|---|
| 1 | GSE120107 | 예 | D0/D4/D7/D15/D30, 각 2반복 | gene count matrix 352 KB | SRA 접근 가능 | 우리 프로토콜과 가장 가까운 주 훈련 기준선; primary-KC 열은 제외 |
| 2 | GSE122383 | 예, hESC H9 | D0/D7/D14/D21/D45, 각 2반복 | processed FPKM 2.4 MB | 12표본 FASTQ 약 150 GB; known-day 10표본만 약 117 GB | 독립 PSC longitudinal validation; D14/D21 공백 보강 |
| 3 | GSE144241 | 예, hESC H1 | D0/D1/D4/D6/D8/D11/D26, 각 1표본 | quantile-normalized matrix 1.6 MB; annotation 7 MB | raw CEL archive 26 MB | 초기 전환점 보조; 반복이 없어 단독 훈련·오차평가 금지 |
| 4 | GSE147206 | 예, hESC/iPSC skin organoid | D6/D29/D48, 5 libraries | 10x processed archive 457 MB; 해제·분석 약 1.7 GB | SRA raw 접근 가능 | scRNA 세포형 특이성, 이질성, off-target 검증; 정확한 2D 날짜 검증은 불가 |
| 5 | GSE287810 | 예, iPSC-KC | D29-D32 endpoint, 6표본 | raw-count matrix 407 KB | GEO/SRA 접근 가능 | 후기 iKC endpoint와 배지 효과 확인; longitudinal 모델에는 부족 |
| 선택 | GSE108248 ATAC-seq | 예 | D0/D7/D14/D21/D43/H9KC 등 | bigWig/peak, GEO archive 6.1 GB | SRA raw 접근 가능 | TFAP2C→TP63 전환의 chromatin 근거; qPCR 날짜 모델의 직접 입력은 아님 |

## 이번에 실제로 확인한 새 자료

### GSE122383

- processed workbook: 27,055행, 16열
- 정상 PSC known-day 열: H9 D0/D7/D14/D21/D45, 각 2반복
- 별도 H9KC 두 열은 선택·확장된 keratinocyte 상태이며 정확한 날짜 anchor로 사용하지 않음
- 15-gene PSC core panel 중 14개가 사용 가능
- 한 replicate trajectory로 다른 replicate를 예측한 2-fold 결과:
  - 평균 절대오차 3일
  - 중앙 절대오차 2.5일
  - 90%가 실제 날짜 ±4일 이내
  - 최대오차 14일: D45 한 표본이 D31로 예측됨
- GSE120107과 비교 가능한 core marker 14개는 모두 시간에 따른 증가/감소 방향이 일치함

해석: PSC→epidermal 단계축은 재현되지만 후기 분화 속도는 protocol/batch 의존적이다. D45를 우리 D35로 그대로 치환하지 않는다.

### GSE144241

- hESC H1 D0/D1/D4/D6/D8/D11/D26
- processed series matrix와 raw CEL이 모두 공개
- D0-D8에는 SB431542, CHIR99021, BMP4, DAPT를 단계적으로 사용하고 D9-D11에는 low-calcium medium과 BMP4/DAPT/EGF를 사용함
- 우리 프로토콜과 약물 조합은 다르지만 PSC pluripotency 소실, surface ectoderm, TP63/KRT5/KRT14 전환 순서를 보조할 수 있음
- 시점당 한 표본이므로 분산이나 오차를 추정할 수 없음

## processed-first 권고

현재 목적에는 raw FASTQ를 먼저 받을 필요가 없다.

1. GSE120107 processed counts로 주 모델을 학습한다.
2. GSE122383 processed FPKM로 완전 독립 PSC 검증을 한다.
3. GSE144241 processed microarray로 D1/D6/D8/D11 초기 전환을 확인한다.
4. GSE147206 processed scRNA-seq에서 epithelial subset과 off-target를 확인한다.
5. 우리 랩 known-day qPCR이 생기면 공개 RNA-seq 값과 직접 합치지 않고, gene-wise 방향을 맞춘 뒤 batch 단위 calibration/test로 나눈다.

raw reprocessing은 다음 경우에만 추천한다.

- GSE120107과 GSE122383을 동일한 genome annotation 및 동일한 quantifier로 재정량해야 할 때
- isoform 수준의 TP63/keratin 분석이 필요할 때
- processed FPKM의 오래된 hg19/TopHat-Cufflinks 처리 영향을 제거해야 할 때

GSE122383 FASTQ 전체를 재처리한다면 다운로드 약 150 GB, 압축 해제·중간파일·정량 결과를 포함한 작업공간 300-500 GB를 잡는 것이 안전하다. Salmon 계열 정량은 일반 워크스테이션에서도 가능하지만, STAR 정렬은 더 많은 RAM과 중간 저장공간이 필요하다.

## 남은 실제 공백

- 우리 프로토콜과 동일한 조건의 D31-D35 PSC 유래 실측값
- D15-D30 사이를 반복 측정한 2D PSC 자료
- 공개 qPCR과 우리 qPCR 사이의 platform calibration

이 공백은 큰 raw 파일을 다시 처리해도 생기지 않는다. 최종적으로는 우리 랩의 D0/D7/D14-15/D21-23/D30/D35 known-day qPCR을 독립 분화 batch 3회 이상 확보해야 해결된다.

## 출처

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120107
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122383
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144241
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147206
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287810
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108248
