# PSC-to-keratinocyte D0-D35 reference

공개 bulk/scRNA-seq 자료와 2D RA/BMP4-EGF-CaCl2 분화 프로토콜을 이용해, 향후 우리 랩 qPCR 표본을 분화 단계와 예상 일자에 매핑하기 위한 reference를 구축한다.

현재 버전은 **실측 공개 시점 사이를 연결한 모델 reference**이다. D0-D35의 모든 날짜가 직접 측정되었다는 뜻은 아니다.

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
| D31-D35 | late maintenance / maturation 방향 | 날짜 신뢰도 낮음 |

## 사용한 공개자료

| 데이터 | 역할 | 핵심 제한 |
|---|---|---|
| GSE120107 | D0, D4, D7, D15, D30 bulk time-course backbone | 각 시점 2반복, D15-D30 간격이 큼 |
| GSE147206 | D6, D29, D48 scRNA-seq 세포형 특이성 | 3D organoid, mixed lineage, D48 한 cell line |
| GSE155816 | primary keratinocyte basal-to-differentiation 보완 | PSC 분화 날짜가 아니라 passage 변화 |
| GSE287810 | D29-D32 iPSC-derived KC endpoint/배지 비교 | endpoint 중심 |
| GSE98483 | primary KC D0-D7 maturation | PSC 분화 날짜와 직접 대응하지 않음 |
| GSE73305 | calcium-induced primary KC maturation | 시점당 생물학적 반복 없음 |

대용량 공개 matrix와 FASTQ는 저장소에 포함하지 않는다. 각 분석 스크립트가 기대하는 accession과 파일명만 기록한다.

## D0-D35 일자별 모델

GSE120107의 D0, D4, D7, D15, D30 평균 발현을 관측 anchor로 사용하고, 그 사이를 유전자별·module별 piecewise-linear interpolation으로 연결했다.

D35에는 같은 자료의 primary keratinocyte를 **상태 anchor**로 배치했다. 이것은 실제 Day 35 측정값이 아니다. D31-D35는 primary-KC-like 방향을 나타내지만 달력 날짜 신뢰도는 낮다.

### 근거 수준

- High: D0, D4, D7, D15, D30 bulk 실측 anchor
- Medium-high: D6, D29 cross-protocol scRNA 지원
- Medium: D5, D8, D14, D23 프로토콜 전환점
- Medium-low: 실측 anchor 사이의 일반 보간일
- Low: D31-D35 후기 상태 보간

## marker panel 초안

- pluripotency: `POU5F1`, `NANOG`
- surface ectoderm/progenitor: `TFAP2A`, `KRT18`, `KRT19`, `TP63`
- basal/epidermal commitment: `KRT5`, `KRT14`, `ITGA6`
- maturation: `KRT1`, `IVL`, `SPRR1B`, `FLG`, `LOR`
- 보조: `DSG1`, `TGM1`, `ABCA12`, `KRT10`
- off-target 감시: `COL1A1`, `PAX6`

GSE147206에서 KRT10은 D6에서도 약 80% 세포에 검출되고 비상피 cluster에도 넓게 나타나 ambient RNA 가능성이 컸다. 따라서 KRT10 단독으로 단계나 세포형을 판정하지 않는다.

반대로 GSE147206에서 거의 검출되지 않았던 IVL과 FLG는 GSE155816에서 충분히 검출되었다. 두 donor 모두 passage 2에서 passage 5로 갈 때 IVL, SPRR1B, FLG, LOR, KRT1이 증가했다. 후기 marker 부재는 marker 자체의 문제라기보다 organoid 자료의 후기 성숙도 또는 library 특성으로 해석한다.

## 데이터 충분성 판단

- processed matrix만으로 marker 방향성, 세포형 특이성, 데이터 공백 평가는 가능했다.
- GSE147206은 QC 후 44,448세포여서 세포 수는 충분했다.
- 그러나 세포 수가 곧 biological replicate 수는 아니다.
- D6 epithelial 비율은 WA25 56.8%, DSP 87.4%로 cell-line 차이가 컸다.
- 정확한 우리 프로토콜 Day 예측에는 known-day qPCR calibration sample이 필요하다.

권장 자체 표본은 D0, D5-D7, D14-D15, D21-D23, D28-D30, D35이며 서로 독립적인 분화 3회 이상이 바람직하다.

## 일자별 모델 자체 검증

15-gene core panel을 사용해 GSE120107에서 reference 날짜를 다시 예측했다.

| 검증 방식 | 평균 절대오차 | ±2일 이내 | 의미 |
|---|---:|---:|---|
| Cross-replicate trajectory | 0.50일 | 100% | 한 replicate로 만든 기준선이 다른 replicate의 관측 anchor를 잘 구분 |
| Leave-one-sample-out | 0.58일 | 100% | 같은 시점의 다른 replicate가 있을 때 안정적 |
| Leave-one-stage-out | 5.13일 | 37.5% | 시점 전체가 없으면 일별 보간 정확도가 크게 감소 |

시점 전체를 숨겼을 때 D4는 D6-D7, D7은 D6, D15는 D20, D30은 D18로 예측됐다. 따라서 현재 모델은 관측 anchor와 가까운 상태를 찾는 데는 유용하지만, D15-D30 사이의 모든 날짜를 정확한 달력 날짜처럼 해석하면 안 된다.

실제 결과는 `best reference day + compatible range + biological stage + fit score` 형식으로 보고한다. 후기 표본은 예를 들어 `best reference day D24; compatible range D20-D29; maturing keratinocyte`처럼 표현한다.

### 훈련셋과 검증셋 분리

시점당 replicate가 2개뿐이므로 무작위 80:20 분할 대신 complete-time-course 2-fold 검증을 사용한다.

- Fold 1: replicate 1의 6개 시점을 training, replicate 2의 6개 시점을 validation
- Fold 2: replicate 2의 6개 시점을 training, replicate 1의 6개 시점을 validation

Fold 1의 평균 절대오차는 0.67일, Fold 2는 0.33일이었다. 두 fold 모두 최대 오차는 2일이고 모든 검증 표본이 실제 anchor의 ±2일 이내였다. 외부 데이터셋은 모델 fitting에서 제외하고 세포형 특이성, endpoint 일치, maturation 방향 검증에만 사용한다.

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
    ├── GSE147206_qpcr_candidate_epithelial_specificity.csv
    ├── GSE155816_marker_passage_effects.csv
    ├── daily_mapping_validation_predictions.csv
    ├── daily_mapping_validation_metrics.csv
    ├── validation_summary_ko.md
    ├── internal_train_validation_split.csv
    ├── internal_cross_replicate_fold_metrics.csv
    ├── external_validation_manifest.csv
    └── train_validation_design_ko.md
```

## 참고자료

- Ali and Abdelalim, Directed differentiation of human pluripotent stem cells into epidermal keratinocyte-like cells, STAR Protocols (2022), https://doi.org/10.1016/j.xpro.2022.101613
- GSE120107: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120107
- GSE147206: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147206
- GSE155816: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155816
- GSE287810: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287810
- GSE98483: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98483
- GSE73305: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73305
