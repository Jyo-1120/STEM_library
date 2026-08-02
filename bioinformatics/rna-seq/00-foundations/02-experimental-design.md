# 실험 설계와 biological replicate

## 1. 한 문장 요약

좋은 RNA-seq 분석은 통계 프로그램이 아니라, 연구 질문에 맞는 조건과 충분한 biological replicate를 정하는 단계에서 시작한다.[^hbc-design]

## 2. 학습 목표

- Biological replicate와 technical replicate를 구분한다.
- Replicate가 조건 내부의 생물학적 변이를 추정하는 데 필요한 이유를 설명한다.
- Confounding과 batch effect의 위험을 설명한다.

## 3. 핵심 개념

### Biological replicate

같은 조건에서 독립적으로 얻은 생물학적 표본이다. 조건 내부의 실제 생물학적 변이를 추정하는 데 사용한다.[^hbc-design]

### Technical replicate

동일한 생물학적 표본을 반복 측정하거나 여러 sequencing lane으로 나눈 경우다. 측정 과정의 변이를 확인할 수 있지만, 독립적인 biological replicate를 대신하지 못한다.[^hbc-design]

### Confounding

관심 조건과 다른 변수가 완전히 겹쳐 효과를 분리할 수 없는 상태다. 예를 들어 control을 모두 첫날 처리하고 treatment를 모두 둘째 날 처리하면 condition과 processing day가 겹친다.

## 4. 실험 전 확인할 질문

- 비교하려는 조건은 무엇인가?
- 표본은 서로 독립적인가?
- 조건마다 biological replicate가 있는가?
- Batch에 모든 조건이 고르게 포함되어 있는가?
- 성별, genotype, 처리 시간처럼 함께 기록해야 할 변수가 있는가?
- Paired design인지 독립 표본 설계인지 명확한가?

## 5. Metadata 예시

```csv
sample_id,condition,batch,subject
control_1,control,1,S01
treated_1,treated,1,S02
control_2,control,2,S03
treated_2,treated,2,S04
```

Batch가 존재한다면 실제 연구 질문과 자료 구조를 검토하여 `design = ~ batch + condition`처럼 모형에 포함할 수 있다. 다만 design formula만으로 이미 완전히 confounded된 실험을 복구할 수는 없다.[^deseq2-design]

## 6. 스스로 확인하기

1. 한 표본을 세 lane에서 sequencing하면 biological replicate가 3개인가?
2. Control과 treatment가 서로 다른 날짜에만 처리되었다면 무엇이 문제인가?
3. Metadata를 sequencing 이후가 아니라 실험 설계 단계에서 준비해야 하는 이유는 무엇인가?

## 7. 참고자료

[^hbc-design]: Harvard Chan Bioinformatics Core, [Experimental design considerations](https://hbctraining.github.io/Intro-to-DGE/lessons/experimental_planning_considerations.html), accessed 2026-08-02.
[^deseq2-design]: Bioconductor, [DESeq2 vignette: Multi-factor designs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs), accessed 2026-08-02.

