# 훈련셋·검증셋 분리 설계

## 내부 검증

GSE120107의 PSC 유래 D0, D4, D7, D15, D30에는 각각 두 replicate가 있다. primary keratinocyte 두 표본은 날짜 모델에서 제외한다. PSC 표본 수가 10개뿐이므로 무작위 80:20 분할 대신, 전체 시간축을 보존하는 2-fold cross-replicate 검증을 사용한다.

### Fold 1

- 훈련: 각 시점 replicate 1, 총 5표본
- 검증: 각 시점 replicate 2, 총 5표본

### Fold 2

- 훈련: 각 시점 replicate 2, 총 5표본
- 검증: 각 시점 replicate 1, 총 5표본

각 fold에서 검증 replicate의 어떤 표본도 reference 생성에 사용하지 않는다.

## 내부 검증 결과

15-gene core panel에서 두 fold 모두 모든 검증 표본이 실제 anchor의 ±2일 이내에 들어왔다.

세부 수치는 `internal_cross_replicate_fold_metrics.csv`에 저장한다. 이 성능은 관측된 anchor 상태를 다른 replicate에서 다시 찾는 능력을 평가하며, 측정되지 않은 중간 날짜의 정확도를 보장하지 않는다.

## 외부 검증셋

훈련에는 GSE120107만 사용하고 다음 데이터셋은 모델 fitting에 사용하지 않는다.

| 데이터셋 | 검증 역할 | 정확한 날짜 검증 가능 여부 |
|---|---|---|
| GSE122383 | hESC D0/D7/D14/D21/D45 독립 longitudinal 검증 | 가능하나 프로토콜 속도 차이를 함께 해석 |
| GSE144241 | hESC D0/D1/D4/D6/D8/D11/D26 초기 전환점 보조 | 제한적: 시점당 1표본 |
| GSE147206 | 세포형 특이성 및 D6/D29 상태 일치 | 불가: 3D organoid·다른 프로토콜 |
| GSE287810 | D29-D32 iKC endpoint 일치 | 제한적: 배지와 프로토콜 차이 |
| GSE155816 | basal-to-differentiation marker 방향 | 불가: primary KC passage 변화 |
| GSE98483 | primary KC 성숙 방향 | 불가: PSC 분화가 아님 |
| GSE73305 | calcium-induced terminal maturation 방향 | 불가: PSC 분화가 아님·시점당 반복 없음 |

외부자료에서는 raw expression 값을 직접 합치지 않고 marker 방향, gene rank, module score를 검증한다.

## 최종 모델 생성

1. 위 두 fold로 panel과 알고리즘을 확정한다.
2. 분석 방법을 고정한 뒤 GSE120107의 PSC 유래 10표본 전체로 public reference를 다시 학습한다.
3. GSE122383을 독립 PSC longitudinal validation으로 사용하고, GSE144241로 초기 marker 전환 방향을 확인한다.
4. 우리 랩 known-day qPCR은 처음에는 완전한 holdout으로 유지한다.
5. 충분한 자체 표본이 모인 뒤에만 일부를 calibration에 포함하고 별도 batch를 최종 test set으로 남긴다.

## 우리 랩 데이터 권장 분리

독립 분화가 3회라면 batch 단위로 나눈다.

- batch 1-2: calibration/training
- batch 3: 완전한 test set

같은 differentiation batch의 서로 다른 well을 train과 test에 나누면 데이터 누출이 생길 수 있으므로 피한다.
