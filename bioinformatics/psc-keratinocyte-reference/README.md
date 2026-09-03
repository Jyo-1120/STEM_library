# PSC-to-keratinocyte D0-D35 reference

공개 bulk/scRNA-seq 자료와 2D RA/BMP4-EGF-CaCl2 분화 프로토콜을 이용해, 향후 우리 랩 qPCR 표본을 분화 단계와 예상 일자에 매핑하기 위한 reference를 구축한다.

현재 버전은 **PSC 유래 실측 공개 시점 사이를 연결한 모델 reference**이다. primary keratinocyte는 날짜 모델에서 제외하며, D0-D35의 모든 날짜가 직접 측정되었다는 뜻은 아니다.

## 목표

우리 랩에서 keratinocyte 분화를 수행하고 qPCR을 측정했을 때 다음을 함께 제시하는 것을 목표로 한다.

1. 가장 가까운 예상 일자
2. 가능한 일자 범위
3. 분화 단계
4. reference와의 일치도
5. off-target 또는 marker 불일치 경고

## 프로토콜 일자 축

| 일자 | 예상 단계 | 주요 처리 |
|---|---|---|
| D0 | pluripotent PSC / 분화 시작 | UCM + RA + BMP4 |
| D1-D4 | ectoderm induction | RA + BMP4 |
| D5 | epithelial induction | 계대, ROCK inhibitor |
| D6-D7 | surface ectoderm 전환 | UCM:N2 전환 후 N2 + RA/BMP4 |
| D8 | progenitor induction 시작 | RA/BMP4 제거, EGF 시작 |
| D9-D13 | early keratinocyte progenitor | N2 + EGF |
| D14 | progenitor checkpoint | p63/KRT18 확인, 계대 |
| D15-D22 | progenitor expansion / basal commitment | N2 + EGF |
| D23-D29 | calcium-responsive maturation | 1.2 mM CaCl2 |
| D30 | mature iKC checkpoint | mature marker 확인 |
| D31-D35 | post-D30 maintenance | 공개 PSC anchor 없음; D30 carry-forward |

## 사용한 공개자료

| 데이터 | 역할 | 핵심 제한 |
|---|---|---|
| GSE120107 | D0, D4, D7, D15, D30 bulk time-course backbone | 각 시점 2반복, D15-D30 간격이 큼 |
| GSE122383 | hESC D0, D7, D14, D21, D45 독립 RNA-seq 시계열 검증 | 후기 분화 속도가 프로토콜 의존적 |
| GSE144241 | hESC D0, D1, D4, D6, D8, D11, D26 microarray | 시점당 1표본; 초기 전환 보조만 가능 |
| GSE147206 | D6, D29, D48 scRNA-seq 세포형 특이성 | 3D organoid, mixed lineage, D48 한 cell line |
| GSE155816 | primary keratinocyte basal-to-differentiation 보완 | PSC 분화 날짜가 아니라 passage 변화 |
| GSE287810 | D29-D32 iPSC-derived KC endpoint/배지 비교 | endpoint 중심 |
| GSE98483 | primary KC D0-D7 maturation | PSC 분화 날짜와 직접 대응하지 않음 |
| GSE73305 | calcium-induced primary KC maturation | 시점당 생물학적 반복 없음 |

대용량 공개 matrix와 FASTQ는 저장소에 포함하지 않는다. 각 분석 스크립트가 기대하는 accession과 파일명만 기록한다.

## D0-D35 일자별 모델

GSE120107의 PSC 유래 D0, D4, D7, D15, D30 평균 발현만 관측 anchor로 사용하고, 그 사이를 유전자별·module별 piecewise-linear interpolation으로 연결했다. 같은 자료의 primary keratinocyte 두 표본은 제외했다.

D31-D35에는 공개 PSC 실측 anchor가 없으므로 D30 값을 그대로 유지한다. 이 구간은 `unanchored`이며 날짜 예측에 사용하지 않는다.

### 근거 수준

- High: D0, D4, D7, D15, D30 bulk 실측 anchor
- Medium-high: D6, D29 cross-protocol scRNA 지원
- Medium: D5, D8, D14, D23 프로토콜 전환점
- Medium-low: 실측 anchor 사이의 일반 보간일
- Low/unanchored: D31-D35 D30 carry-forward

## marker panel 초안

- pluripotency: `POU5F1`, `NANOG`
- surface ectoderm/progenitor: `TFAP2A`, `KRT18`, `KRT19`, `TP63`
- basal/epidermal commitment: `KRT5`, `KRT14`, `ITGA6`
- maturation: `KRT1`, `IVL`, `SPRR1B`, `ABCA12`, `DSG1`, `TGM1`
- terminal maturation QC only: `FLG`, `LOR`
- 보조: `KRT10`
- off-target 감시: `COL1A1`, `PAX6`

GSE147206에서 KRT10은 D6에서도 약 80% 세포에 검출되고 비상피 cluster에도 넓게 나타나 ambient RNA 가능성이 컸다. 따라서 KRT10 단독으로 단계나 세포형을 판정하지 않는다.

반대로 GSE147206에서 거의 검출되지 않았던 IVL과 FLG는 GSE155816에서 충분히 검출되었다. 두 donor 모두 passage 2에서 passage 5로 갈 때 IVL, SPRR1B, FLG, LOR, KRT1이 증가했다. 후기 marker 부재는 marker 자체의 문제라기보다 organoid 자료의 후기 성숙도 또는 library 특성으로 해석한다.

## 데이터 충분성 판단

- processed matrix만으로 marker 방향성, 세포형 특이성, 데이터 공백 평가는 가능했다.
- GSE122383 processed FPKM는 2.4 MB이며, raw FASTQ는 전체 약 150 GB이다.
- GSE144241 processed matrix는 1.6 MB, raw CEL archive는 26 MB이다.
- GSE147206은 QC 후 44,448세포여서 세포 수는 충분했다.
- 그러나 세포 수가 곧 biological replicate 수는 아니다.
- D6 epithelial 비율은 WA25 56.8%, DSP 87.4%로 cell-line 차이가 컸다.
- 정확한 우리 프로토콜 Day 예측에는 known-day qPCR calibration sample이 필요하다.

권장 자체 표본은 D0, D5-D7, D14-D15, D21-D23, D28-D30, D35이며 서로 독립적인 분화 3회 이상이 바람직하다.

## 훈련·검증 구조

- 주 훈련: GSE120107 PSC D0/D4/D7/D15/D30, 10표본
- 내부 검증: replicate trajectory 2-fold, 평균 절대오차 1.1일, 전 표본 ±2일 이내
- 독립 PSC 검증: GSE122383 H9 D0/D7/D14/D21/D45, 10표본
- GSE122383 내부 replicate 검증: 평균 절대오차 3일, 90%가 ±4일 이내
- 초기 시점 보조: GSE144241 D0/D1/D4/D6/D8/D11/D26

GSE120107과 GSE122383에서 비교 가능한 core marker 14개는 모두 시간에 따른 증가/감소 방향이 일치했다. 다만 D45 한 표본이 D31로 예측되어 후기 날짜는 프로토콜 간 직접 이전하지 않는다.

## qPCR에 적용할 때

RNA-seq log2CPM과 qPCR delta-Ct를 직접 비교하지 않는다. 우리 랩 known-day qPCR 자료가 확보되면 다음 순서로 calibration한다.

1. reference gene을 이용해 delta-Ct 계산
2. 발현 방향을 맞추기 위해 `-delta-Ct` 사용
3. 유전자별 중심과 scale 표준화
4. 일자별 reference와 거리 또는 상관계수 계산
5. 최적 일자와 함께 허용 범위, 단계명, 일치도 출력

## 저장소 구조

```text
bioinformatics/psc-keratinocyte-reference/
├── README.md
├── scripts/
│   ├── analyze_gse120107.R
│   ├── analyze_gse122383.R
│   ├── analyze_gse144241.R
│   ├── analyze_gse147206_scRNA.R
│   ├── analyze_gse155816_marker_detection.py
│   ├── analyze_maturation_references.R
│   ├── build_daily_reference.R
│   ├── validate_daily_reference.R
│   └── define_train_validation_sets.R
└── results/
    ├── daily_stage_reference.csv
    ├── model_anchors.csv
    ├── GSE120107_qpcr_candidate_summary.csv
    ├── GSE122383_internal_cross_replicate_metrics.csv
    ├── GSE122383_cross_study_core_gene_direction.csv
    ├── GSE144241_core_gene_trends.csv
    ├── psc_only_dataset_accessibility_ko.md
    ├── GSE147206_qpcr_candidate_epithelial_specificity.csv
    └── GSE155816_marker_passage_effects.csv
```

## 참고자료

- Ali and Abdelalim, Directed differentiation of human pluripotent stem cells into epidermal keratinocyte-like cells, STAR Protocols (2022), https://doi.org/10.1016/j.xpro.2022.101613
- GSE120107: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120107
- GSE122383: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122383
- GSE144241: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144241
- GSE147206: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147206
- GSE155816: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155816
- GSE287810: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287810
- GSE98483: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98483
- GSE73305: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73305
