## RNA-seq count matrix와 metadata 불러오기

## 파일 경로는 실제 프로젝트 구조에 맞게 수정합니다.
counts <- read.csv(
  "data/counts.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata.csv",
  row.names = 1
)

## 기본 무결성 확인
stopifnot(!anyNA(counts))
stopifnot(all(counts >= 0))
stopifnot(setequal(colnames(counts), rownames(metadata)))

## Metadata를 count matrix 열 순서에 맞춥니다.
metadata <- metadata[colnames(counts), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))

