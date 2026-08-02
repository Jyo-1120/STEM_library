# Count normalization

## 1. 한 문장 요약

Count normalization은 sequencing depth와 RNA composition 같은 관심 밖의 기술적 차이를 조정하여 sample 사이의 count를 비교할 수 있게 만드는 과정이다.[^hbc-normalization]

## 2. 학습 목표

- Library size와 RNA composition의 영향을 설명한다.
- DESeq2 median-of-ratios 방법의 큰 흐름을 설명한다.
- Normalization과 VST/rlog 변환을 구분한다.

## 3. 왜 필요한가?

관측 count는 실제 RNA 발현뿐 아니라 sequencing depth와 sample composition의 영향을 받는다. 따라서 raw count의 크기만 직접 비교하면 기술적 차이를 생물학적 차이로 잘못 해석할 수 있다.[^hbc-normalization]

## 4. DESeq2의 기본 접근

DESeq2는 sample별 size factor를 추정하며, 기본적으로 median-of-ratios 방법을 사용한다.[^deseq2-size-factor]

개념적으로는 다음 순서로 이해할 수 있다.

1. Gene마다 sample 전체의 대표값을 계산한다.
2. 각 count를 해당 gene의 대표값과 비교한다.
3. Sample별 ratio들의 중앙값을 size factor로 사용한다.
4. Raw count를 size factor로 조정해 normalized count를 얻는다.

## 5. R 코드

```r
library(DESeq2)

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
```

전체 차등발현 분석에서 `DESeq(dds)`를 실행하면 필요한 size factor 추정 단계도 수행된다.[^hbc-deseq2]

## 6. Normalization과 transformation

| 작업 | 목적 | 대표 사용처 |
|---|---|---|
| Size-factor normalization | Sample 간 기술적 규모 차이 보정 | DESeq2 모형, normalized count 확인 |
| VST/rlog | 평균에 따른 분산 변화를 완화 | PCA, clustering, heatmap |

VST나 rlog로 변환한 값을 일반적인 DESeq2 차등발현 검정의 입력으로 다시 넣지 않는다.[^deseq2-transform]

## 7. 흔한 오해

- “정규화하면 모든 sample의 총합이 같아진다”라고 단순화하지 않는다.
- Normalization은 batch effect 제거와 같은 의미가 아니다.
- Normalized count와 VST/rlog 값은 같은 데이터가 아니다.
- TPM으로 변환하는 것이 모든 sample 간 차이를 해결하는 것은 아니다.

## 8. 확인 문제

1. Sequencing depth가 큰 sample은 항상 생물학적으로 발현량이 높은가?
2. DESeq2가 sample별 size factor를 계산하는 이유는 무엇인가?
3. PCA에는 raw count보다 VST/rlog 값을 사용하는 이유가 무엇인가?

## 9. 참고자료

[^hbc-normalization]: Harvard Chan Bioinformatics Core, [Count normalization with DESeq2](https://hbctraining.github.io/Intro-to-DGE/lessons/02_DGE_count_normalization.html), accessed 2026-08-02.
[^deseq2-size-factor]: Bioconductor, [DESeq2 vignette: The DESeq2 model](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseq2-model), accessed 2026-08-02.
[^hbc-deseq2]: Harvard Chan Bioinformatics Core, [Description of steps for DESeq2](https://hbctraining.github.io/Intro-to-DGE/lessons/04b_DGE_DESeq2_analysis.html), accessed 2026-08-02.
[^deseq2-transform]: Bioconductor, [DESeq2 vignette: Variance stabilizing transformation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation), accessed 2026-08-02.

