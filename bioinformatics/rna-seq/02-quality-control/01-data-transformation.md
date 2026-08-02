# 데이터 변환: VST와 rlog

## 1. 한 문장 요약

VST와 rlog는 count 데이터에서 평균이 클수록 분산도 커지는 경향을 완화하여 sample 간 거리, PCA, clustering 같은 탐색적 분석을 돕는다.[^deseq2-transform]

## 2. 왜 필요한가?

Raw count에서는 발현량이 높은 gene의 분산이 큰 경향이 있어 거리 계산과 PCA가 일부 고발현 gene에 지나치게 영향을 받을 수 있다. DESeq2는 이를 위한 변환으로 VST와 rlog를 제공한다.[^deseq2-transform]

## 3. R 코드

```r
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)
```

일반적으로 큰 데이터셋에서는 VST가 rlog보다 빠르다. `blind` 설정은 변환에서 실험 설계를 얼마나 고려할지 결정하므로, 목적과 데이터 구조에 맞춰 공식 vignette의 설명을 확인해야 한다.[^deseq2-blind]

## 4. 중요한 구분

VST/rlog는 sample QC와 시각화를 위한 변환이다. 차등발현 검정 자체는 변환된 값이 아니라 count에 대한 음이항 모형을 사용한다.[^deseq2-transform]

## 5. 확인 문제

1. VST는 sequencing depth 보정을 의미하는가?
2. VST 값을 `DESeqDataSetFromMatrix()`의 count 입력으로 사용해도 되는가?
3. `blind = TRUE`와 `blind = FALSE`는 어떤 상황에서 검토해야 하는가?

## 6. 참고자료

[^deseq2-transform]: Bioconductor, [DESeq2 vignette: Variance stabilizing transformation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation), accessed 2026-08-02.
[^deseq2-blind]: Bioconductor, [DESeq2 vignette: Blind dispersion estimation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#blind-dispersion-estimation), accessed 2026-08-02.

