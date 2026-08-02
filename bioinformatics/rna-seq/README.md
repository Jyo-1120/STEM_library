# Bulk RNA-seq와 차등발현 분석

Harvard Chan Bioinformatics Core(HBC)의 **Introduction to Differential Gene Expression Analysis** 과정을 뼈대로 삼아, bulk RNA-seq의 gene-level 차등발현 분석을 한국어로 다시 정리합니다.[^hbc-course]

이 노트는 강의 내용을 번역해서 옮기는 것이 아니라 다음 질문에 스스로 답할 수 있도록 만드는 것을 목표로 합니다.

- 이 분석 단계는 왜 필요한가?
- 입력 데이터와 출력 데이터는 무엇인가?
- 그래프와 통계 결과에서 무엇을 판단해야 하는가?
- 실제 연구 데이터에 적용할 때 무엇을 주의해야 하는가?

## 학습 범위

```text
FASTQ
  ↓
품질 확인 및 정량화
  ↓
gene × sample count matrix  ← 이 노트의 주된 시작점
  ↓
DESeq2 객체 생성 및 정규화
  ↓
샘플 수준 품질 확인
  ↓
통계 모형과 비교 조건 설정
  ↓
차등발현 분석
  ↓
결과 시각화와 기능 분석
```

HBC 과정은 Salmon으로 정량화한 결과를 `tximport`로 가져오는 예제를 사용하며, 이후 DESeq2를 이용한 gene-level 분석을 설명합니다.[^hbc-setup][^tximport]

## 학습 순서

### 0. 기초

1. [RNA-seq 전체 흐름](00-foundations/01-rnaseq-overview.md)
2. [실험 설계와 biological replicate](00-foundations/02-experimental-design.md)
3. [FASTQ에서 count matrix까지](00-foundations/03-from-fastq-to-counts.md)

### 1. Count 데이터

1. [Count matrix 이해하기](01-count-data/01-count-matrix.md)
2. RNA-seq count 분포
3. [Count normalization](01-count-data/03-normalization.md)

### 2. 샘플 수준 품질 확인

1. [데이터 변환: VST와 rlog](02-quality-control/01-data-transformation.md)
2. [PCA](02-quality-control/02-pca.md)
3. Sample correlation
4. Hierarchical clustering

### 3. 차등발현 분석

1. Design formula
2. DESeq2 workflow
3. Wald test와 contrast
4. 다중검정 보정
5. LRT와 time-course 분석

### 4. 시각화

1. MA plot
2. Volcano plot
3. Heatmap

### 5. 기능 분석

1. Gene annotation
2. Over-representation analysis
3. GSEA

## 학습 체크리스트

- [ ] RNA-seq 전체 분석 흐름을 설명할 수 있다.
- [ ] Raw count, normalized count, TPM의 용도를 구분할 수 있다.
- [ ] Biological replicate가 필요한 이유를 설명할 수 있다.
- [ ] DESeq2 normalization의 목적을 설명할 수 있다.
- [ ] VST/rlog와 normalization을 구분할 수 있다.
- [ ] PCA와 sample correlation 결과를 해석할 수 있다.
- [ ] 연구 질문에 맞는 design formula와 contrast를 설정할 수 있다.
- [ ] `log2FoldChange`, `pvalue`, `padj`를 해석할 수 있다.
- [ ] Wald test와 LRT의 사용 목적을 구분할 수 있다.
- [ ] ORA와 GSEA의 차이를 설명할 수 있다.

## 출처 표시 원칙

- 본문에서 근거가 필요한 설명 뒤에 각주를 표시합니다.
- 직접 인용은 따옴표와 출처를 함께 표시합니다.
- 원문을 요약한 경우에도 원자료 링크를 남깁니다.
- 강의 자료와 공식 소프트웨어 문서를 우선 사용합니다.
- 블로그나 질의응답은 보조 자료로만 사용합니다.

## 주요 참고자료

[^hbc-course]: Harvard Chan Bioinformatics Core, [Introduction to Differential Gene Expression Analysis: Lessons](https://hbctraining.github.io/Intro-to-DGE/schedule/links-to-lessons.html), accessed 2026-08-02.
[^hbc-setup]: Harvard Chan Bioinformatics Core, [Set up and overview for gene-level differential expression analysis](https://hbctraining.github.io/Intro-to-DGE/lessons/01b_DGE_setup_and_overview.html), accessed 2026-08-02.
[^tximport]: Bioconductor, [Importing transcript abundance with tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html), accessed 2026-08-02.

- [DESeq2 공식 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [HBC Intro-to-DGE Quarto 버전](https://hbctraining.github.io/Intro-to-DGE-Quarto/)
- [HBC 과정 전체 일정과 단원 설명](https://hbctraining.github.io/Intro-to-DGE/schedule/)
- [참고자료 모음](references/useful-links.md)

## 라이선스 및 인용 안내

HBC 교육자료는 CC BY 4.0 조건으로 공개되어 있습니다. 이 저장소에서 해당 자료를 바탕으로 작성한 문서에는 원저자와 원문 링크를 표시합니다. 자세한 저자 및 인용 정보는 [HBC 과정 소개](https://hbctraining.github.io/Intro-to-DGE-Quarto/)에서 확인할 수 있습니다.

