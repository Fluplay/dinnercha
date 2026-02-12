function [A,B] = dimensionForce(A,id1,B,id2)
    % <Description>

    % make A's id1 leg and B's id2 leg contractable, and fill in random data

    % <Input>

    % A : [QSpace tensor]
    % B : [QSpace tensor]
    % id1 id2 : [numeric] leg number to force same dimension

    % <Output>
    % A, B : [QSpace tensor]

    % Written by Subin Kim (01.02.2025)
    QA = A.Q{id1};
    QB = B.Q{id2};
    [Qcom, indices] = compareVectors(QA, QB);
    for i = 1:size(Qcom,1)
        dimA = size(A.data{indices(i,1)},id1);
        dimB = size(B.data{indices(i,2)},id2);
        [mDim, index] = min([dimA dimB]); % minimum dimension of common quantum number
        if index == 1
            sectors = findVectorIndices(QB, Qcom(i,:)); % Qcom(i,L) is common qunatum number
            for j = sectors
                B.data{j} = 0.1+ rand(size(sliceTensorDimension(B.data{j},id2,mDim)));
            end
        else
            sectors = findVectorIndices(QA, Qcom(i,:)); % Qcom(i,L) is common qunatum number
            for j = sectors
                A.data{j} = 0.1 + rand(size(sliceTensorDimension(A.data{j},id1,mDim)));
            end
        end
    end
end


function [result, indices] = compareVectors(QA, QB)
    [nA, mA] = size(QA);
    [nB, mB] = size(QB);
    
    if mA ~= mB
        error('QA와 QB의 열 수가 같아야 합니다.');
    end
    
    result = [];
    indices = [];
    
    for i = 1:nA
        for j = 1:nB
            if isequal(QA(i,:), QB(j,:))
                result = [result; QA(i,:)];
                indices = [indices; i, j];
                break;  % 매치를 찾았으므로 내부 루프를 종료
            end
        end
    end
    
    if isempty(result)
        result = zeros(0, mA);
    end
end

function indices = findVectorIndices(A, targetVector)
    % A: nxm double 배열
    % targetVector: 길이 m의 벡터
    
    % A와 targetVector의 열 수가 같은지 확인
    [n, m] = size(A);
    if length(targetVector) ~= m
        error('targetVector의 길이가 A의 열 수와 같아야 합니다.');
    end
    
    % 결과를 저장할 변수 초기화
    indices = [];
    
    % A의 각 행을 targetVector와 비교
    for i = 1:n
        if isequal(A(i,:), targetVector)
            indices = [indices i];
        end
    end
    
    % 결과가 없을 경우 빈 배열 반환
    if isempty(indices)
        indices = [];
    end
end

function slicedTensor = sliceTensorDimension(tensor, id, mDim)
    % tensor: 입력 텐서
    % id: 슬라이스할 차원의 인덱스
    % mDim: 슬라이스할 범위 (1:mDim)

    % 텐서의 차원 수 확인
    dims = ndims(tensor);
    
    % 동적 인덱싱을 위한 셀 배열 생성
    idx = repmat({':'}, 1, dims);
    
    % 지정된 차원에 대해 인덱스 범위 설정
    idx{id} = 1:mDim;
    
    % 텐서 슬라이싱
    slicedTensor = tensor(idx{:});
end