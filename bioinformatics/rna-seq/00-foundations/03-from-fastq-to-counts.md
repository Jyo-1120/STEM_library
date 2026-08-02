# FASTQ에서 count matrix까지

## 1. 한 문장 요약

차등발현 분석 전에는 sequencing read의 품질을 확인하고, reference에 대한 정렬 또는 정량화를 통해 gene별 count를 생성해야 한다.[^hbc-workflow]

## 2. 두 가지 대표 경로

### Alignment 기반

```text
FASTQ → STAR/HISAT2 → BAM → featureCounts → gene count matrix
```

### Transcript 정량화 기반

```text
FASTQ → Salmon/Kallisto → transcript abundance → tximport → gene-level matrix
```

`tximport`는 transcript-level abundance, estimated counts, transcript length를 가져와 downstream gene-level 분석에 사용할 행렬로 요약한다.[^tximport]

## 3. DESeq2에 전달하기 전에 확인할 것

- Reference genome과 annotation 버전
- Gene ID 체계: Ensembl ID, Entrez ID, gene symbol
- Count matrix의 sample 이름과 metadata의 sample 이름
- Raw gene count인지, Salmon 등의 estimated count인지
- TPM/FPKM과 count를 혼동하지 않았는지

## 4. 범위 메모

이 저장소의 첫 학습 단계에서는 HBC의 Salmon 예제 결과에서 시작한다. FASTQ 품질 관리와 정량화는 이후 별도 실습으로 확장한다.

## 5. 참고자료

[^hbc-workflow]: Harvard Chan Bioinformatics Core, [RNA-seq workflow](https://hbctraining.github.io/Intro-to-DGE/lessons/01a_RNAseq_processing_workflow.html), accessed 2026-08-02.
[^tximport]: Bioconductor, [Importing transcript abundance with tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html), accessed 2026-08-02.

