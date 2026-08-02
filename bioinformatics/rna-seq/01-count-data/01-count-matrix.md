# Count matrix 이해하기

## 1. 한 문장 요약

Count matrix는 각 gene에서 관측 또는 추정된 read 수를 sample별로 정리한 행렬이며, DESeq2의 통계 모형과 직접 연결되는 핵심 입력이다.[^deseq2-input]

## 2. 기본 구조

| gene_id | control_1 | control_2 | treated_1 | treated_2 |
|---|---:|---:|---:|---:|
| GeneA | 105 | 98 | 240 | 221 |
| GeneB | 0 | 2 | 1 | 0 |
| GeneC | 1200 | 1310 | 1180 | 1255 |

- 각 행은 동일한 annotation에 정의된 gene을 나타내야 한다.
- 각 열은 하나의 sample과 대응해야 한다.
- 일반적인 count matrix 입력에서는 음수와 결측값을 허용하지 않는다.
- 정렬 후 직접 계산한 count를 사용할 때 DESeq2는 정수형 count를 기대한다.[^deseq2-input]

## 3. 서로 구분해야 할 값

| 값 | 주된 용도 |
|---|---|
| Raw count | Count 기반 차등발현 모형의 입력 |
| Normalized count | Sample 간 발현 패턴 탐색과 일부 시각화 |
| VST/rlog 값 | PCA, clustering 등 sample QC |
| TPM | Sample 내부에서 상대적인 abundance 확인 등에 사용 |

DESeq2의 차등발현 검정에는 미리 정규화한 값이 아니라 원 count 또는 지원되는 transcript quantification 결과를 사용한다. DESeq2가 모형 내부에서 size factor를 추정하기 때문이다.[^deseq2-why-counts]

## 4. 데이터 무결성 확인 예시

```r
stopifnot(!anyNA(counts))
stopifnot(all(counts >= 0))
stopifnot(identical(colnames(counts), rownames(metadata)))
```

## 5. 스스로 확인하기

1. Count matrix의 열 이름과 metadata의 행 이름이 달라도 되는가?
2. TPM 테이블을 DESeq2의 일반적인 raw count 입력처럼 사용하면 안 되는 이유는 무엇인가?
3. Gene ID가 중복된다면 어떤 원인을 먼저 확인해야 하는가?

## 6. 참고자료

[^deseq2-input]: Bioconductor, [DESeq2 vignette: Count matrix input](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input), accessed 2026-08-02.
[^deseq2-why-counts]: Bioconductor, [DESeq2 vignette: Why un-normalized counts?](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts), accessed 2026-08-02.

