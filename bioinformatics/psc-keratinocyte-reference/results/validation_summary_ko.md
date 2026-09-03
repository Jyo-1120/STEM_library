# D0-D35 일자별 reference 자체 검증

## 결론

현재 15-gene core panel은 **이미 실측 anchor가 있는 상태를 다시 찾아내는 능력은 좋다.** 그러나 특정 시점 전체를 숨겼을 때 후기 날짜를 복원하는 능력은 낮아, D15-D30 사이의 일별 값은 정확한 달력 날짜가 아니라 `reference day index`로 해석해야 한다.

## 사용한 core panel

`POU5F1`, `NANOG`, `TFAP2A`, `KRT18`, `KRT19`, `TP63`, `KRT5`, `KRT14`, `ITGA6`, `KRT1`, `IVL`, `SPRR1B`, `FLG`, `LOR`, `ABCA12`

KRT10은 GSE147206에서 ambient RNA 가능성이 있어 제외했다. DSG1, TGM1, COL1A1, PAX6를 추가한 expanded panel도 비교했으나 core panel보다 성능이 낮았다.

## 검증 1: cross-replicate trajectory

한 replicate의 D0/D4/D7/D15/D30/primary-KC 전체로 기준선을 만들고, 다른 replicate의 모든 시점을 예측했다. 그 반대 방향도 반복했다.

- 평균 절대오차: 0.5일
- 중앙 절대오차: 0일
- 모든 표본이 실제 anchor의 ±2일 이내
- 평균 최적 상관계수: 0.947
- marker bootstrap 범위가 모든 실제 anchor를 포함

이 결과는 core panel이 공개자료의 알려진 상태를 안정적으로 구분한다는 뜻이다. 다만 모든 anchor 시점이 학습 reference에 있으므로, 측정되지 않은 중간 날짜의 정확도를 직접 증명하지는 않는다.

## 검증 2: leave-one-sample-out

각 표본을 하나씩 제외하고 나머지 11개 표본으로 기준선을 다시 만든 뒤 제외한 표본의 날짜를 예측했다.

- 평균 절대오차: 0.58일
- 중앙 절대오차: 0일
- 모든 표본이 ±2일 이내
- 평균 최적 상관계수: 0.949

이 결과도 replicate 간 재현성이 좋다는 근거이지만, 해당 시점의 다른 replicate가 학습에 남아 있다는 제한이 있다.

## 검증 3: leave-one-stage-out interpolation

D4, D7, D15 또는 D30의 두 replicate를 모두 제거하고, 남은 앞뒤 시점만으로 숨긴 시점을 예측했다.

- 평균 절대오차: 5.13일
- 중앙 절대오차: 4일
- ±2일 이내: 37.5%
- ±4일 이내: 50%
- 평균 최적 상관계수: 0.833

| 숨긴 실제 시점 | 예측 결과 |
|---|---|
| D4 | D6-D7 |
| D7 | D6 |
| D15 | D20 |
| D30 | D18 |

D4-D7 구간은 비교적 가까이 복원됐지만, D15와 D30 사이의 변화는 직선 보간으로 복원되지 않았다. 특히 D30을 숨기면 D18로 예측되어 후기 분화가 비선형임을 보여준다.

## 실제 qPCR에 적용할 출력 형식

현재 모델은 다음처럼 출력하는 것이 적절하다.

- `best reference day`: 가장 가까운 일자 template
- `day range`: marker bootstrap 및 근접 template 범위
- `stage`: pluripotent, surface ectoderm, progenitor, basal-like, maturing KC
- `fit`: reference와의 상관계수 또는 거리
- `warning`: off-target, marker discordance, low-confidence late interpolation

예를 들어 후기 표본을 `D24`로 단정하기보다 `best reference day D24, compatible range D20-D29, maturing keratinocyte`처럼 보고한다.

## 현재 해석 규칙

- D0-D7: 날짜 수준의 참고값으로 사용 가능하나 범위를 함께 제시
- D8-D15: 단계 판정은 가능, 정확한 날짜는 중간 신뢰도
- D16-D29: basal-to-maturation 상태축으로 사용, 날짜 범위는 넓게 제시
- D30: 직접 anchor가 있으므로 D30-like 상태 판정 가능
- D31-D35: primary-KC-like 성숙 방향이며 실제 달력 날짜 판정으로 사용하지 않음

## 다음 보강점

우리 프로토콜에서 D21-D23 또는 D25-D28 known-day qPCR 표본을 하나만 추가해도 D15-D30 사이의 가장 큰 공백을 나눌 수 있다. 최소 설계는 D0, D7, D14/15, D23, D30, D35이며, 독립 분화 3회 이상을 권장한다.
