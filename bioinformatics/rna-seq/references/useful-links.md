# RNA-seq 참고자료

가능하면 교육자료보다 공식 package 문서를 최종 기준으로 확인합니다. Package는 계속 업데이트될 수 있으므로 분석 기록에는 R 및 package 버전을 함께 남깁니다.

## 기반 교육과정

- Harvard Chan Bioinformatics Core, [Introduction to Differential Gene Expression Analysis](https://hbctraining.github.io/Intro-to-DGE/schedule/links-to-lessons.html)
- Harvard Chan Bioinformatics Core, [Intro-to-DGE Quarto version](https://hbctraining.github.io/Intro-to-DGE-Quarto/)
- 원본 GitHub 저장소, [hbctraining/Intro-to-DGE](https://github.com/hbctraining/Intro-to-DGE)

## 공식 package 문서

- Bioconductor, [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- Bioconductor, [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- Bioconductor, [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
- Bioconductor, [tximport vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
- Bioconductor, [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

## 논문

- Love MI, Huber W, Anders S. [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8). Genome Biology. 2014.
- Soneson C, Love MI, Robinson MD. [Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](https://doi.org/10.12688/f1000research.7563.2). F1000Research. 2015.

## 출처 작성 예시

```markdown
DESeq2는 count 데이터에 음이항 일반화 선형모형을 적용한다.[^deseq2]

[^deseq2]: Love MI, Huber W, Anders S.
    [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8).
    Genome Biology. 2014.
```

