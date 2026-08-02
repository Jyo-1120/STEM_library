# RNA-seq 전체 흐름

## 1. 한 문장 요약

RNA-seq 분석은 RNA에서 얻은 sequencing read를 유전자 또는 transcript의 발현량으로 정량화하고, 조건 사이의 차이가 생물학적 변화인지 통계적으로 평가하는 과정이다.[^hbc-workflow]

## 2. 학습 목표

- FASTQ부터 차등발현 결과까지의 흐름을 설명한다.
- 원시 데이터 처리와 count matrix 이후 분석을 구분한다.
- 각 단계에서 만들어지는 대표적인 파일을 설명한다.

## 3. 전체 과정

| 단계 | 대표 입력 | 주요 작업 | 대표 출력 |
|---|---|---|---|
| 실험 설계 | 연구 질문 | 조건, 반복, batch 결정 | Sample sheet |
| Sequencing | RNA library | 염기서열 판독 | FASTQ |
| Raw QC | FASTQ | Read 품질과 adapter 확인 | QC report |
| 정량화 | FASTQ, reference | Alignment 또는 lightweight mapping | Transcript/gene counts |
| 데이터 구성 | Count 파일, metadata | 샘플과 조건 연결 | Count matrix |
| Sample QC | Count matrix | PCA, correlation, clustering | QC plots |
| 차등발현 분석 | Count matrix, metadata | 정규화, 모형 적합, 검정 | DEG results |
| 기능 분석 | DEG 또는 gene ranking | ORA, GSEA | Pathway/GO results |

HBC 과정은 raw read로부터 count를 생성하는 전처리 흐름을 소개하지만, 실습의 중심은 Salmon 정량화 결과를 가져온 이후의 gene-level 차등발현 분석이다.[^hbc-workflow][^hbc-setup]

## 4. 가장 중요한 두 데이터

### Count matrix

- 행: gene
- 열: sample
- 값: 해당 sample에서 관측되거나 추정된 gene별 count

### Sample metadata

- 각 sample의 조건과 실험 정보를 기록한다.
- Count matrix의 열 이름과 metadata의 sample 이름이 정확히 대응해야 한다.
- Condition뿐 아니라 batch, sex, genotype, time point 등 분석에 필요한 변수를 포함할 수 있다.

## 5. 스스로 확인하기

1. FASTQ 파일을 DESeq2에 직접 입력하지 않는 이유는 무엇인가?
2. Count matrix와 metadata는 어떤 기준으로 연결되는가?
3. Sample QC는 차등발현 검정 전후 중 언제 수행해야 하는가?

## 6. 참고자료

[^hbc-workflow]: Harvard Chan Bioinformatics Core, [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://hbctraining.github.io/Intro-to-DGE/lessons/01a_RNAseq_processing_workflow.html), accessed 2026-08-02.
[^hbc-setup]: Harvard Chan Bioinformatics Core, [Set up and overview for gene-level differential expression analysis](https://hbctraining.github.io/Intro-to-DGE/lessons/01b_DGE_setup_and_overview.html), accessed 2026-08-02.

