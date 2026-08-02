# PCA로 sample 관계 확인하기

## 1. 한 문장 요약

PCA는 많은 gene의 발현 정보를 소수의 축으로 요약하여 sample 사이의 주요 변이 구조와 잠재적인 이상치를 탐색하는 방법이다.[^hbc-pca]

## 2. PCA에서 확인할 것

- 같은 조건의 biological replicate가 가까운가?
- 관심 조건에 따라 sample이 분리되는가?
- Batch, sex, sequencing date 같은 다른 변수가 분리를 설명하는가?
- 다른 sample들과 매우 멀리 떨어진 sample이 있는가?

PCA에서 분리되어 보인다는 사실만으로 sample을 바로 제외해서는 안 된다. Metadata, library QC, read quality, mapping/quantification 결과와 함께 원인을 조사해야 한다.[^hbc-sample-qc]

## 3. R 코드

```r
vsd <- vst(dds, blind = TRUE)

plotPCA(
  vsd,
  intgroup = c("condition", "batch")
)
```

직접 꾸미기 위해 PCA 데이터를 추출할 수도 있다.

```r
pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "batch"),
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))
```

## 4. 해석할 때 주의할 점

- PC1과 PC2가 전체 변이를 모두 보여주는 것은 아니다.
- 축의 방향과 부호 자체에는 절대적인 생물학적 의미가 없다.
- PCA는 차등발현 유전자를 판정하는 통계 검정이 아니다.
- 조건별 분리가 작다고 해서 차등발현 유전자가 반드시 없는 것은 아니다.

## 5. 기록할 내용

PCA 그림만 저장하지 말고 다음을 함께 기록한다.

```text
- 사용한 데이터 변환:
- 표시한 metadata 변수:
- PC1 설명 분산:
- PC2 설명 분산:
- 조건별 clustering 여부:
- Batch pattern 여부:
- 이상치 후보와 추가 확인 사항:
```

## 6. 확인 문제

1. PCA에서 같은 조건의 sample 하나가 멀리 떨어졌다면 어떤 정보를 추가로 확인해야 하는가?
2. Sample이 PC1에서 batch별로 나뉜다면 design formula에서 무엇을 검토해야 하는가?
3. PC1과 PC2에서 조건별 분리가 없으면 DEG 분석을 중단해야 하는가?

## 7. 참고자료

[^hbc-pca]: Harvard Chan Bioinformatics Core, [Principal Component Analysis](https://hbctraining.github.io/Intro-to-DGE/lessons/principal_component_analysis.html), accessed 2026-08-02.
[^hbc-sample-qc]: Harvard Chan Bioinformatics Core, [Sample-level QC](https://hbctraining.github.io/Intro-to-DGE/lessons/03_DGE_QC_analysis.html), accessed 2026-08-02.

