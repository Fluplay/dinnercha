function O = getRandTensor(A, legA, B, legB, D)
    if dim(A,legA)*dim(B,legB) < D
        O = getIdentity(A,legA,B,legB);
        D3 = dim(O,3);
        for i = 1:numel(O.data)
            bond3 = size(O.data{i},3);
            bond3Aft = ceil(bond3/D3*(D));
            bond3Aft = min([bond3Aft bond3]);
            tensorSize = [size(O.data{i},1) size(O.data{i},2) bond3Aft];
            O.data{i} = randn(tensorSize);
        end
        return
    end

    % get random tensor of shape [A1 x B1 x D], with direction in-in-out
    
    % 작은 크기의 A와 B 생성 (차원은 1로, 구조는 유지)
    A_small = reduceToDim1(A);
    B_small = reduceToDim1(B);
    
    % 작은 텐서들을 사용하여 결합된 구조 생성
    AB_small = getIdentity(A_small, legA, B_small, legB, [1 2 3]);
    
    % 최종 결과 텐서의 구조 준비
    O = AB_small;
    
    % 3번째 leg에 존재하는 고유한 대칭 섹터 식별
    [uniqueQLabels, sectorDims] = identifyUniqueThirdLegSectors(AB_small, A, legA, B, legB, D);
    
    % 각 블록마다 적절한 크기의 랜덤 데이터 생성
    for i = 1:numel(O.data)
        % 첫 번째, 두 번째 leg의 실제 차원 추출
        dimA = getActualLegDim(A, legA, O.Q{1}(i,:));
        dimB = getActualLegDim(B, legB, O.Q{2}(i,:));
        
        % 현재 블록의 3번째 leg 대칭 라벨
        qLabel3rd = O.Q{3}(i,:);
        
        % 해당 대칭 섹터의 차원 찾기
        for j = 1:size(uniqueQLabels, 1)
            if all(qLabel3rd == uniqueQLabels(j,:))
                dimD = sectorDims(j);
                break;
            end
        end
        
        % 랜덤 데이터 생성
        O.data{i} = randn(dimA, dimB, dimD);
    end
    
    % 적절한 itag 설정
    O = setitags(O, {A.info.itags{legA}, B.info.itags{legB}, 'rand*'});
end

% 3번째 leg의 고유 대칭 섹터를 식별하고 각 섹터에 차원 할당
function [uniqueQLabels, sectorDims] = identifyUniqueThirdLegSectors(AB_small, A, legA, B, legB, totalDim)
    % 3번째 leg의 모든 대칭 라벨 수집
    allQLabels3rd = AB_small.Q{3};
    
    % 고유한 대칭 라벨 찾기
    [uniqueQLabels, ~, ic] = unique(allQLabels3rd, 'rows');
    numUniqueSectors = size(uniqueQLabels, 1);
    
    % 각 고유 섹터에 대한 평균 차원 비율 계산
    sectorEstimatedDims = zeros(numUniqueSectors, 1);
    sectorBlockCounts = zeros(numUniqueSectors, 1);
    
    for i = 1:numel(AB_small.data)
        % 현재 블록의 대칭 라벨
        qLabelsA = AB_small.Q{1}(i,:);
        qLabelsB = AB_small.Q{2}(i,:);
        sectorIdx = ic(i);
        
        % 현재 블록에 대한 3번째 leg 차원 추정
        dimEstimate = estimateThirdLegDim(qLabelsA, qLabelsB, uniqueQLabels(sectorIdx,:), A, legA, B, legB);
        
        % 해당 섹터의 추정 차원에 추가
        sectorEstimatedDims(sectorIdx) = sectorEstimatedDims(sectorIdx) + dimEstimate;
        sectorBlockCounts(sectorIdx) = sectorBlockCounts(sectorIdx) + 1;
    end
    
    % 각 섹터의 평균 추정 차원 계산
    for i = 1:numUniqueSectors
        if sectorBlockCounts(i) > 0
            sectorEstimatedDims(i) = sectorEstimatedDims(i) / sectorBlockCounts(i);
        else
            sectorEstimatedDims(i) = 1; % 기본값
        end
    end
    
    % 섹터 차원 비율 계산
    sectorRatios = sectorEstimatedDims / sum(sectorEstimatedDims);
    
    % 전체 차원 D를 섹터 비율에 따라 분배
    sectorDims = max(1, round(totalDim * sectorRatios));
    
    % 총합이 정확히 D가 되도록 조정
    while sum(sectorDims) ~= totalDim
        if sum(sectorDims) < totalDim
            [~, idx] = max(sectorRatios .* (1 - sectorDims/sum(sectorDims)));
            sectorDims(idx) = sectorDims(idx) + 1;
        else
            [~, idx] = max(sectorDims .* (1 - sectorRatios));
            sectorDims(idx) = sectorDims(idx) - 1;
        end
    end
end

% 3번째 leg 차원 추정 함수
function dimEstimate = estimateThirdLegDim(qLabelsA, qLabelsB, qLabelsAB, A, legA, B, legB)
    % 대칭 유형 확인
    qtype = A.info.qtype;
    isSU2xSU2 = contains(qtype, 'SU2') && length(strfind(qtype, 'SU2')) >= 2;
    
    if isSU2xSU2
        % SU(2)×SU(2) 대칭성의 경우
        % A와 B의 다중도 차원 계산
        dimA_multiplet = calculateSU2Dimension(qLabelsA);
        dimB_multiplet = calculateSU2Dimension(qLabelsB);
        
        % 3번째 leg 다중도 차원 계산
        dimAB_multiplet = calculateSU2Dimension(qLabelsAB);
        
        % 해당 블록의 실제 차원
        dimA_actual = getActualLegDim(A, legA, qLabelsA);
        dimB_actual = getActualLegDim(B, legB, qLabelsB);
        
        % A와 B의 실제 차원 비율
        ratioA = dimA_actual / dimA_multiplet;
        ratioB = dimB_actual / dimB_multiplet;
        
        % 3번째 leg의 추정 차원
        dimEstimate = dimAB_multiplet * sqrt(ratioA * ratioB);
    else
        % 기타 대칭성의 경우 (간단한 방식)
        dimEstimate = 1;
    end
end

% 이전에 정의된 보조 함수들
function X_small = reduceToDim1(X)
    X_small = X;
    for i = 1:numel(X.data)
        sz = size(X.data{i});
        new_sz = ones(size(sz));
        new_sz(1:length(sz)) = 1;
        X_small.data{i} = ones(new_sz);
    end
end

function dim = getActualLegDim(X, leg, qLabels)
    for i = 1:size(X.Q{leg}, 1)
        if all(X.Q{leg}(i,:) == qLabels)
            dim = size(X.data{i}, leg);
            return;
        end
    end
    dim = 1;
end

function dim = calculateSU2Dimension(qLabels)
    dim = 1;
    for i = 1:length(qLabels)
        dim = dim * (qLabels(i) + 1);
    end
end