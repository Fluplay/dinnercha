function idx = circularIndex(n, maxN)
    % CIRCULARINDEX 1부터 maxN까지의 순환 인덱스를 반환합니다
    %   n: 현재 인덱스
    %   maxN: 최대 인덱스 값
    %   반환값: 1과 maxN 사이의 순환 인덱스
    
    % 입력 유효성 검사
    validateattributes(n, {'numeric'}, {'integer'});
    validateattributes(maxN, {'numeric'}, {'integer', 'positive'});
    
    % mod 함수를 사용한 순환 인덱스 계산
    idx = mod(n - 1, maxN) + 1;
end